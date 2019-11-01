%% cosmosFS - cosmos Flight Simulator
%% emulates: 
%% -realtime code execution through timed execution and pauses
%% -execution on different satellites through parallel programming

%% to do:
%% -physics: SRP, J2, aerodynamics
%% -run cases to prove extended PL through solar pressure
%% -count memory and flops, estimate storage requirements, conclude whether possible to run on Arduino
%% -HIL
%% -elliptical orbit
%% -start further workers to emulate OBC, GPS for mode switching
%% -parellel-run orbitDetermination
%% -clarify what data of the header goes to each satellite
%% -use of Matlab timer function, maybe
%% -one-RW control law
%% -respect P/L attitude needs

%% Revision history:
%% 29/09/2019:  Ivanov case works
%% 26/09/2019:  cycle works, computational algorithm incomplete

close all;clear all;clc;%#ok<CLALL>
oldpath = path; path(oldpath,'..\matlabfunctions\')
oldpath = path; path(oldpath,'..\cosmosSupport\')

%% symbolic names of initial conditions and desired statevector functions
%sstInitialFunction=@IRSRendezvousInitial;
%sstInitialFunction=@IvanovFormationFlightInitial;
sstInitialFunction=@cluxterInitial;

%% actual initial conditions of ODE, altitude is not used here therefore the ~
[~,ns,~,panels,rho,v,radiusOfEarth,~,mu,satelliteMass,panelSurface,...
  sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite]=sstInitialFunction(); 

DQ = parallel.pool.DataQueue;
afterEach(DQ,@disp);
parpool(ns);

startTime=posixtime(datetime('now')); %% posixtime, i.e. seconds
accelerationFactor=10000;
maxOrbits=120;

%% initial idx and altitude
idx=120;
altitude=340000;  

%% data that will later be per satellite and therefore inside SPMD loop
orbitSection    =2;         %degree
orbitSections   =[1:orbitSection:360];
orbitCounter    =0;
error           =zeros(6,ns);
sst             =zeros(9,1);
sstDesired      =zeros(6,1);
sstOld          =zeros(9,1);
refPosChangeTemp=zeros(3,1);

%% these variables are for plotting. They will not be used for operational software
SST_PP            =zeros(9,1);
REFPOSCHANGE_PP   =zeros(3,1);
TIME_PP           =0;


%% non-gravitational perturbations
wind          =windOn*rho/2*v^2*[-1 0 0]';
sunlight      =sunOn*2*4.5e-6*[0 0 -1]' ;       %% at reference location, needs to be rotated later

%refsurf=panelSurface*panels(1);
refSurf=panelSurface*panels(3);

%% forcevector determination and its angular granulaty
alphas            =0:deltaAngle:360;     %% roll
betas             =0:deltaAngle:180;     %% pitch 
gammas            =0:deltaAngle:360;     %% yaw

aeropressureforcevector  =aeropressureforcevectorfunction(wind,panelSurface,panels(1),...
                                                          panels(2),panels(3),alphas,betas,gammas);
solarpressureforcevector =solarpressureforcevectorfunction(sunlight,panelSurface,panels(1),...
                                                          panels(2),panels(3),alphas,betas,gammas);

spmd(ns) %% create satellite instances
%%   this is code for each satellite, some of the code of above will come down here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  alive=1;
  send(DQ,strcat('No. ',num2str(labindex),' is alive.'))
  while alive
    %% satellites are alive but not doing anything
    [goFoFli,batteryOK]=getMode(maxOrbits,orbitCounter,DQ);
    if batteryOK
      %switch on GPS
    end
    SST_PP   =sst;
    
    while (goFoFli==1 || goFoFli==2)  %% orbit loop
      orbitCounter=orbitCounter+1;
      startOrbit=now; %% posixtime, i.e. seconds
      %send(DQ,strcat(num2str(labindex),': orbittimer begin: ',num2str(posixtime(datetime('now'))-startTime)));

      [meanMotion,meanAnomalyFromAN,altitude,sst]=whereInWhatOrbit(sst,altitude,(idx-1)/size(orbitSections,2));
      %% settings for control algorithm
      [P,IR,A,B]=riccatiequation(meanMotion/180*pi);  
      %% 
      %send(DQ,strcat(num2str(labindex),': MeanMotion: ',num2str(meanMotion)));
      %send(DQ,strcat(num2str(labindex),': meanAnomalyFromAN: ',num2str(meanAnomalyFromAN)));
      
      %% wait until end of orbit section
      idx = find(orbitSections >= meanAnomalyFromAN,1,'first');
      idx=idx+1;

      pause((orbitSections(idx)-meanAnomalyFromAN)/meanMotion/accelerationFactor);

      while (goFoFli==1 || goFoFli==2) && idx<=size(orbitSections,2) %% orbit sections loop
        startSection=now; %% determine cycle start time in order allow to subtract cycle duration from waiting
        %% set attitude computed in last iteration
        setAttitude();
        %% compute attitude for next section
        %% determine desired trajectory
        sstDesired=sstDesiredFunction(orbitSections(idx)/meanMotion,meanMotion/180*pi,labindex,goFoFli);
        %[sstDesired]=cluxterDesired(timetemptemp,MeanMotion,i)
        %% determine error
        error(1:6,labindex)=sst(1:6)-sstDesired(1:6);
        if not(masterSatellite)
          
          %! ISL error
          for i=1:ns
            if i~=labindex
              labSend(error(:,labindex),i); %% send my error to i
              error(:,i) = labReceive(i); %% receive i's error, error contains the error vector for each i
            end
          end
          averageError=zeros(6,ns);
          %! compute average error
          for i=1:ns %% compute average error
            averageError(:,i)=averageError(:,i)+error(:,i)/ns; %% averageError contains columns of average error computed on each i
          end
          %! ISL averageerror
          for i=1:ns
            if i~=labindex
              labSend(averageError(:,i),i);
              averageError(:,i)= labReceive(i);
            end
          end
          %% check whether all errors are the same, this error check is relatively memory intensive, worth it?
          errorSum=0;
          for i=2:ns
            for j=1:6
              errorSum=errorSum+averageError(j,1)-averageError(j,i);
            end
          end
          if errorSum~=0
            send(DQ,strcat(num2str(labindex),': error sum NOK:',num2str(errorSum) ));
          end
          %! assign average error, i.e. compute final error
          error(2:6,labindex) =error(2:6,labindex)-averageError(2:6,labindex);
          error(1,labindex)   =error(1,labindex)-max(error(1,:));
        end
        
        %% compute attitude
        sstOld=sst;
        [sst,~]=HCWEquation(...
          IR,P,A,B,orbitSection/meanMotion,sstOld,error(:,labindex),...
          sqrt(wind(1)^2+wind(2)^2+wind(3)^2),sqrt(sunlight(1)^2+sunlight(2)^2+sunlight(3)^2),...
          alphas,betas,gammas,aeropressureforcevector,solarpressureforcevector,sstOld(7),...
          sstOld(8),sstOld(9),refSurf,satelliteMass,wakeAerodynamics,masterSatellite,...
          orbitSections(idx)/meanMotion,radiusOfEarth,altitude,meanMotion/180*pi);
        sst=sst';
        TIME_PP=[TIME_PP TIME_PP(end)+orbitSection/meanMotion];  %% for plotting
        
        if not(masterSatellite)
          %% send reference position to all non 1-satellites and receive by them          
          refPosChange=zeros(3,1);
          if labindex==1
            refPosChange(1:3)=sst(1:3)-sstOld(1:3);
          end
          if labindex==1
            for i=2:ns
              tag=1000000*i+10000*idx+100*orbitCounter+1;
              labSend(refPosChange,i,tag);
            end
          end
          for i=2:ns
            if labindex==i
              tag=1000000*i+10000*idx+100*orbitCounter+1;            
              refPosChange= labReceive(1,tag);
            end
          end                 
          
          %! move coordinate system
          sst(1:3)=sst(1:3)-refPosChange(1:3);
          REFPOSCHANGE_PP=[REFPOSCHANGE_PP refPosChange]; %% for plotting
        end
        %% set counter and go condition
        idx=idx+1;
        [goFoFli,batteryOK]=getMode(maxOrbits,orbitCounter,DQ);      
        %% wait until next section starts
        %send(DQ,strcat(num2str(labindex),': time4section',num2str(orbitSection/meanMotion),' now-start ',num2str(now-start)));
        if idx>3/4*size(orbitSections,2) && batteryOK %% orbit determination will be switch on when in sun light and battery is ok
          orbitDetermination();
        end
        pause(orbitSection/meanMotion/accelerationFactor-(now-startSection));
        SST_PP=[SST_PP sst]; %% for plotting
      end %% orbit sections while loop
      send(DQ,strcat(num2str(labindex),': - no of orbit: ',num2str(orbitCounter),' duration: ',num2str(now-startOrbit),' s -'));  
    end %% orbit while loop 
    send(DQ,strcat('No.',num2str(labindex),' is dead.'))
    if goFoFli~=1 %% if orbit is broken, also break alive loop, this will change later with other conditions
      alive=0;
    end
  end %% alive loop
end %% parallel computation

%% post processing for/and vizualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sstPlot=zeros(size(SST_PP{1},1),ns,size(SST_PP{1},2));
for iii=1:ns
  sstPlot(:,iii,:)=SST_PP{iii};
end

%% results' plotting
u             =zeros(3,ns,size(SST_PP{1},2));
e             =zeros(6,ns,size(SST_PP{1},2));
plotting(sstPlot(7:9,:,:),sstPlot(1:6,:,:),REFPOSCHANGE_PP{1},TIME_PP{1},ns,meanMotion{1}/180*pi,u,e);

delete(gcp)

function setAttitude()
  pause(0.001);
end

function [meanMotion,meanAnomalyFromAN,altitude,sst]=whereInWhatOrbit(sst,altitude,endOfSectionsCycle)
%% this function usually provides a meanAnomalyFromAN=0 and the related meanMotion for a circular orbit
  
  %% info that will come from GPS or TLE
  [GPSTLEdataAvailable,altitudeGPSTLE,meanAnomalyFromANGPSTLE,time]=getGNSSOrTLEdata();
    
  if GPSTLEdataAvailable
    altitude=altitudeGPSTLE;
    meanAnomalyFromAN=meanAnomalyFromANGPSTLE;
    %% compute SST if possible
  else
    if endOfSectionsCycle
      meanAnomalyFromAN=0.01;
    else
      meanAnomalyFromAN=120;
    end    
  end
  %% use sst, meanAnomalyFromAN and altitude either from GPS or from input parameters, %! define rule
  [rho,v,radiusOfEarth,mu,meanMotion]=orbitalproperties(altitude);
  meanMotion=meanMotion/pi*180;
end


function orbitDetermination()
%% switches orbit determination based on GPS or TLE on
%% this function needs to run in parallel to the main thread and must provide trajectory to getGNSSOrTLEdata()
%% how the heck does this work?
  GPSMethod=1;
  if GPSMethod==1
    ;
    %switchGPS(1);
    %getGPSData();
    %computeTrajectoryGPS();
  else
    ;
    %getTLE();
    %computeTrajectoryTLE();
  end
end

function [GPSTLEdataAvailable,altitudeGPSTLE,meanAnomalyFromANGPSTLE,time]=getGNSSOrTLEdata()
%% compute meanAnomalyFromANGPS,altitudeGPS at now from available past GPS or TLE/SGP4 data
  GPSTLEdataAvailable=0;
  altitudeGPSTLE=1;
  meanAnomalyFromANGPSTLE=1;
  time=1;
end

function [goFoFli,batteryOK]=getMode(maxOrbits,orbitCounter,DQ)
  goFoFli=1;
  %% readfile
  readModeFromFile=10; %% default: go to alive loop/mode
  if orbitCounter>=maxOrbits
    send(DQ,strcat(num2str(labindex),': leaving loop - maximum number of orbits reached'))
    goFoFli=10;
  elseif readModeFromFile==0
    send(DQ,strcat(num2str(labindex),': leaving loop - filestop'))
    goFoFli=20;
  else
    if orbitCounter<10
      goFoFli=1;
    else
      goFoFli=2;
    end
  end
  batteryOK=1;
end

function go=setMode(startTime,endTime,DQ)
  go=1;
end

function satellite(ns,i,DQ)
%timeline=0:1:100; 
%error=zeros(ns,1);
%averageError=0;
%  for t=1:size(timeline,2)
%  end
  OK = 9876;
  data = 1234;
  if labindex==1
    send(DQ,labindex);    
    labSend(200,2);
    labSend(300,3);
  elseif labindex==2
    send(DQ,labindex);
    data=labReceive(1);
    send(DQ,data);    
  elseif labindex==3
    send(DQ,labindex);
    data=labReceive(1);
    send(DQ,data);    
  else
    send(DQ,1234567890)
  end 
end


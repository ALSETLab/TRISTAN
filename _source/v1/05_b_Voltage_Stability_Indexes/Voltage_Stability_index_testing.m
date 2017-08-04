clc; clear all;close all
warning off

global simTime faultime deltime Lb samples

% simTime  -  simulation time of the data
% faultime -  Fault occuring time in the system
% deltime  -  Load variation interval
% Lb       -  Load bus
%samples   -  Size of samples fromt the interval data

% Chose between one of the two case studies to understand the  voltage stability function

flag=1;   % 1 for 2-bus test system  
          % 2 for KTH-Nordic 32 bus system for tripping of line 32
          % 3 for KTH-Nordic 32 bus system for tripping of line 42
          
%%          
if flag==1
load 2bus-system;
simTime  = 20;  % Simulation time
faultime = 5;   % Fault occuring time
deltime  = 5;   % Load Variation time interval
Lb=2;           % Load bus to be analyzed
samples=10;     % Size of the sampled interval from data
Limits=[4 4];   % Limits=[Plim,Vlim] in (%),  Percent of Power from the Nose
else if flag==2
load KTH_Nordic_trip_L32;        
simTime  = 40;   % Simulation time
faultime = 10;   % Fault occuring time
deltime  = 10;   % Load Variation time interval
Lb=40;           % Load bus to be analyzed
samples=10;     % Size of the sampled interval from data
Limits=[1 1];    % Limits=[Plim,Vlim] in (%),  Percent of Power from the Nose
else if flag==3    
load KTH_Nordic_trip_L42;        
simTime  = 40;   % Simulation time
faultime = 10;   % Fault occuring time
deltime  = 10;   % Load Variation time interval
Lb=47;           % Load bus to be analyzed
samples=10;     % Size of the sampled interval from data
Limits=[1 1];    % Limits=[Plim,Vlim] in (%),  Percent of Power from the Nose
    end
    end
end
t=ts;
P=-Pinj(:,Lb);  %The Power Injection for the load bus is given as negative in PSAT
Q=-Qinj(:,Lb);  %The Power Injection for the load bus is given as negative in PSAT
V=Vmag(:,Lb);
figure(1000); plot(t,V,'--b')
title('Voltage at bus where the load is connected')

figure(2000); plot(t,P,'--r')
title('Active Power at line where the load is connected')

%% Function for Post-contingency PV curve estimation

global alpha beta   % Alpha and Beta are active and reactive power relation of the load i.e.Q=alpha+Beta*P
[PPostmax,VPostmax,P_pvpost,V_pvpost,state,E_post,X_post,PPost,VPost,alp,bet]=dynamic_voltagestability_postcontingency(t,P,V,Q,faultime,deltime,simTime,samples);

%% Pre-contingency-No load condition
for m=1:length(Lb)
P0(1,m)=0;
Q0(1,m)=0;
if flag==1
V0(1,m)=1.0002;      % Voltage of the load bus in 2 bus system 
else if flag==2
V0(1,m)=0.95856;     % Voltage of the load bus in KTH-Nordic 32 bus system 
else if flag==3
V0(1,m)=0.9979;      % Voltage of the load bus in KTH-Nordic 32 bus system 
    end
    end
end
    
disp('For No load condition the following P Q V are considered:');
disp(P0(1,m));disp(Q0(1,m));disp(V0(1,m));
prompt='Do you wish to change the voltage (Y/N) ? ';
str = input(prompt,'s');
if str=='Y'||str=='y'
prompt='Enter the no load voltage at the load bus: ';
V0(1,m)=input(prompt);  
else
V0(1,m)=V0(1,m);  
end

%% Pre-contingency-fixing the range of maximum active power loading (P_max) and minimum Voltage (V_min)
%For 2-bus test system
if flag==1
 P2(:,m)=[1.5;2.0];
 V2(:,m)=[0.86081;0.78337];
else if flag==2
%For KTH-Nordic 32 bus system
P2(:,m)=[7.5;8.25];
V2(:,m)=[0.87147;0.86286];
    else if flag==3
  P2(:,m)=[9.0;10.0];
  V2(:,m)=[0.89843;0.88784];
        end  
    end
end
P_start(1,m)=P2(1,m);P_end(1,m)=P2(2,m);
V_start(1,m)=V2(1,m);V_end(1,m)=V2(2,m);

disp('The Range of P_max is:');disp(P2(:,m)); 
prompt='Do you wish to change the range of maximum line loading P_max (Y/N) ? ';
str = input(prompt,'s');

if str=='Y'||str=='y'
prompt='Enter the starting value in the range of P_max: ';
P_start(1,m)=input(prompt);
prompt='Enter the ending value in the range of P_max: ';
P_end(1,m)=input(prompt);
end

disp('The Range of V_min considered: ');disp(V2(:,m)); 
prompt='Do you wish to change the range of V_min (Y/N) ?';
str = input(prompt,'s');
if str=='Y'||str=='y'
prompt='Enter the minimum value in the range of V_min';
V_start(:,m)=input(prompt);
prompt='Enter the maximum value in the range of V_min';
V_end(:,m)=input(prompt);
end

P2(:,m)=[P_start(:,m);P_end(:,m)];
V2(:,m)=[V_start(:,m);V_end(:,m)];
Q2(:,m)=alp(m)+bet(m).*P2(:,m);
end
%% Function for Pre-contingency PV curve estimation

[PPremax1,VPremax1,PPremax2,VPremax2,P_pvpre1,V_pvpre1,P_pvpre2,V_pvpre2,states1,states2,E_pre,X_pre,P1,V1]=dynamic_voltagestability_precontingency(t,P,V,Q,P0,Q0,V0,P2,V2,Q2,faultime,alp,bet,samples);

%% Assigning the percentage of limits
for aa=1:length(Lb)
Pmax0_pre(aa)=Limits(1);
Pmax0_post(aa)=Limits(1);
Vmax0_pre(aa)=Limits(2);
end

PPremax=[PPremax1;PPremax2];
VPremax=[VPremax1;VPremax2];
%% Function for caculation of Voltage stability indexes/indices
[SBI1,ABI1,GBI1,SBI2,ABI2,GBI2,Pre1,Post1,Pre2,Post2]=indices_calculation(PPremax,VPremax,Pmax0_pre,Pmax0_post,Vmax0_pre,P0,P1,P2,V0,V1,V2,P_pvpost,V_pvpost,PPost,VPost,PPostmax,VPostmax,P_pvpre1,V_pvpre1,P_pvpre2,V_pvpre2);

%% Displaying the ouputs
clc
disp('The Post-contingency values of the given network');
disp('X_Th :');disp(X_post);
disp('E :');disp(E_post);
disp('alpha:');disp(alpha);
disp('beta:');disp(beta);

disp('The Pre-contingency values of the given network');
disp('X_Th :');disp(X_pre);
disp('E :');disp(E_pre);
disp('The following V_min and P_max values are considered');
disp('Max P Loading  Min.V ');
disp([P2 V2]);

disp('Considering the first pre-contingency PV curve (Lower limit of maximum loading)');
SBI1
ABI1
GBI1

disp('Considering the second pre-contingency PV curve (Upper limit of maximum loading)');
SBI2
ABI2
GBI2

% Outputs
% SBI1,SBI2 - Single Bus Index. Is a matrix and provides:  the distance in pre and post contingency 
%       for each loading level to  the maximum loadability (Pmax) and
%       voltage (Vlim) w.r.t pre-contingency curve 1 & 2 for the given load bus data.
% ABI1,ABI2 - All Bus Index. Is a vector that provides the minimum distance among all buses for each 
%       loading level  to the loadability (Pmax) and the voltage (Vlim) limits w.r.t pre-contingency curve 1 & 2 for the given load bus data.
%       The ABI index helps to identify if a loading level is violating a limit in pre or post contingency.
% GBI1,GBI2 - Is a 2 element vector that provides the overall minimum distance to  the loadability (Pmax)
%       and the over all voltage (Vmax) limits respect to all buses  w.r.t pre-contingency curve 1 & 2 for the given load bus data. 
%       The GBI index indicates if a limit has been violated.
% Pre1 Pre2 Post1 & Post2 - Structure containing the following information. 
%     .limits - row 1 Power limits, row 2 Voltage limits
%     .mean   - row 1, Power mean values for each loading level [Pa,Pb,Pb], where a=low, b=OK, c=limit levels 
%               row 2, Voltage mean values for each loading level [Va,Vb,Vc],where a=low, b=OK, c=limit levels 
%     .deltas - same information as SBI
%     .curve  - column 1 estimated Power, useful to reproduce nose curve
%               column 2 estimated Voltage, useful tor efproduce nose curve
%                    i.e. precontingency PV plot = plot(Pre.curve(:,1),Pre.curve(:,2)) 

return

%% Optional, PLotting Resulting Curves

qq=length(Lb); 

for j=1:qq,
pcent1=20/100;
pcent2=5/100;

Pm1=Pre1(:,:,j).limits(1)+(Pre1(:,:,j).limits(1)*(Limits(1)/100));
xmin1= Pre1(:,:,j).mean(1,1)-( Pre1(:,:,j).mean(1,1)*pcent1);
%xmax= Pre(:,:,j).limits(1)+(Pre(:,:,j).limits(1)*pcent2);
xmax1= Pm1+(Pm1*pcent2);
ymin1= Pre1(:,:,j).limits(2)-(Pre1(:,:,j).limits(2)*pcent1);
ymax1= Pre1(:,:,j).mean(2,1)+(Pre1(:,:,j).mean(2,1)*pcent2);


figure; 
plot(Pre1(:,:,j).mean(1,:),Pre1(:,:,j).mean(2,:),'g*')
title('pre-contingency curve-1')
hold on;plot(Pre1(:,:,j).curve(:,1),Pre1(:,:,j).curve(:,2),'b-')
axis([xmin1,xmax1,ymin1,ymax1])

hline(Pre1(:,:,j).limits(2),'r:',['Vlim1= ',num2str(Pre1(:,:,j).limits(2))]);
vline(Pre1(:,:,j).limits(1),'r:',['Plim1= ',num2str(Pre1(:,:,j).limits(1))]);
vline(Pre1(:,:,j).mean(1,1),'b:',['Pa1= ',num2str(Pre1(:,:,j).mean(1,1))]);
vline(Pre1(:,:,j).mean(1,2),'b:',['Pb1= ',num2str(Pre1(:,:,j).mean(1,2))]);
vline(Pre1(:,:,j).mean(1,3),'b:',['Pc1= ',num2str(Pre1(:,:,j).mean(1,3))]);
hline(Pre1(:,:,j).mean(2,1),'b:',['Va1= ',num2str(Pre1(:,:,j).mean(2,1))]);
hline(Pre1(:,:,j).mean(2,2),'b:',['Vb1= ',num2str(Pre1(:,:,j).mean(2,2))]);
hline(Pre1(:,:,j).mean(2,3),'b:',['Vc1= ',num2str(Pre1(:,:,j).mean(2,3))]);

figure; 
Pm2=Pre2(:,:,j).limits(1)+(Pre2(:,:,j).limits(1)*(Limits(1)/100));
xmin2= Pre2(:,:,j).mean(1,1)-( Pre2(:,:,j).mean(1,1)*pcent1);
%xmax= Pre(:,:,j).limits(1)+(Pre(:,:,j).limits(1)*pcent2);
xmax2= Pm2+(Pm2*pcent2);
ymin2= Pre2(:,:,j).limits(2)-(Pre2(:,:,j).limits(2)*pcent1);
ymax2= Pre2(:,:,j).mean(2,1)+(Pre2(:,:,j).mean(2,1)*pcent2);
 
plot(Pre2(:,:,j).mean(1,:),Pre2(:,:,j).mean(2,:),'b*')
title('pre-contingency curve-2');
hold on;plot(Pre2(:,:,j).curve(:,1),Pre2(:,:,j).curve(:,2),'g-')
axis([xmin2,xmax2,ymin2,ymax2])

hline(Pre2(:,:,j).limits(2),'r:',['Vlim2= ',num2str(Pre2(:,:,j).limits(2))]);
vline(Pre2(:,:,j).limits(1),'r:',['Plim2= ',num2str(Pre2(:,:,j).limits(1))]);
vline(Pre2(:,:,j).mean(1,1),'b:',['Pa2= ',num2str(Pre2(:,:,j).mean(1,1))]);
vline(Pre2(:,:,j).mean(1,2),'b:',['Pb2= ',num2str(Pre2(:,:,j).mean(1,2))]);
vline(Pre2(:,:,j).mean(1,3),'b:',['Pc2= ',num2str(Pre2(:,:,j).mean(1,3))]);
hline(Pre2(:,:,j).mean(2,1),'b:',['Va2= ',num2str(Pre2(:,:,j).mean(2,1))]);
hline(Pre2(:,:,j).mean(2,2),'b:',['Vb2= ',num2str(Pre2(:,:,j).mean(2,2))]);
hline(Pre2(:,:,j).mean(2,3),'b:',['Vc2= ',num2str(Pre2(:,:,j).mean(2,3))]);
 
figure;
plot(Pre1(:,:,j).mean(1,:),Pre1(:,:,j).mean(2,:),'k*')
title('post-contingency with Pre-contingency curve-1')
hold on;plot(Pre1(:,:,j).curve(:,1),Pre1(:,:,j).curve(:,2),'b-')
hold on;plot(Post1(:,:,j).curve(:,1),Post1(:,:,j).curve(:,2),'r-')
hold on;plot(Post1(:,:,j).mean(1,:),Post1(:,:,j).mean(2,:),'c*')
axis([xmin1,xmax1,ymin1,ymax1])

hline(Pre1(:,:,j).limits(2),'b:',['Vlim_pre1= '  ,num2str(Pre1(:,:,j).limits(2))]);
vline(Pre1(:,:,j).limits(1),'b:',['Plim_pre1= ' ,num2str(Pre1(:,:,j).limits(1))]);
vline(Post1(:,:,j).limits(1),'r:',['Plim_post= ',num2str(Post1(:,:,j).limits(1))]);
vline(Pre1(:,:,j).mean(1,1),'k:',['Pa= ',num2str(Post1(:,:,j).mean(1,1))]);
vline(Pre1(:,:,j).mean(1,2),'k:',['Pb= ',num2str(Post1(:,:,j).mean(1,2))]);
vline(Pre1(:,:,j).mean(1,3),'k:',['Pc= ',num2str(Post1(:,:,j).mean(1,3))]);
hline(Post1(:,:,j).mean(2,1),'k:',['Va= ',num2str(Post1(:,:,j).mean(2,1))]);
hline(Post1(:,:,j).mean(2,2),'k:',['Vb= ',num2str(Post1(:,:,j).mean(2,2))]);
hline(Post1(:,:,j).mean(2,3),'k:',['Vc= ',num2str(Post1(:,:,j).mean(2,3))]);

figure;
plot(Pre2(:,:,j).mean(1,:),Pre2(:,:,j).mean(2,:),'k*')
title('post-contingency with Pre-contingency curve-2')

hold on;plot(Pre2(:,:,j).curve(:,1),Pre2(:,:,j).curve(:,2),'g-')
hold on;plot(Post2(:,:,j).curve(:,1),Post2(:,:,j).curve(:,2),'r-')
hold on;plot(Post2(:,:,j).mean(1,:),Post2(:,:,j).mean(2,:),'c*')
axis([xmin2,xmax2,ymin2,ymax2])


hline(Pre2(:,:,j).limits(2),'g:',['Vlim_pre2= '  ,num2str(Pre2(:,:,j).limits(2))]);
vline(Pre2(:,:,j).limits(1),'g:',['Plim_pre2= ' ,num2str(Pre2(:,:,j).limits(1))]);
vline(Post2(:,:,j).limits(1),'r:',['Plim_post= ',num2str(Post2(:,:,j).limits(1))]);
vline(Pre2(:,:,j).mean(1,1),'k:',['Pa= ',num2str(Post2(:,:,j).mean(1,1))]);
vline(Pre2(:,:,j).mean(1,2),'k:',['Pb= ',num2str(Post2(:,:,j).mean(1,2))]);
vline(Pre2(:,:,j).mean(1,3),'k:',['Pc= ',num2str(Post2(:,:,j).mean(1,3))]);
hline(Post2(:,:,j).mean(2,1),'k:',['Va= ',num2str(Post2(:,:,j).mean(2,1))]);
hline(Post2(:,:,j).mean(2,2),'k:',['Vb= ',num2str(Post2(:,:,j).mean(2,2))]);
hline(Post2(:,:,j).mean(2,3),'k:',['Vc= ',num2str(Post2(:,:,j).mean(2,3))]);

figure;
plot(Pre1(:,:,j).mean(1,:),Pre1(:,:,j).mean(2,:),'b*')
hold on;plot(Pre2(:,:,j).mean(1,:),Pre2(:,:,j).mean(2,:),'g*')
title('All the PV curves')
hold on;plot(Pre1(:,:,j).curve(:,1),Pre1(:,:,j).curve(:,2),'b-')
hold on;plot(Pre2(:,:,j).curve(:,1),Pre2(:,:,j).curve(:,2),'g-')
hold on;plot(Post2(:,:,j).curve(:,1),Post2(:,:,j).curve(:,2),'r-')
hold on;plot(Post1(:,:,j).mean(1,:),Post1(:,:,j).mean(2,:),'b*')
hold on;plot(Post2(:,:,j).mean(1,:),Post2(:,:,j).mean(2,:),'g*')
axis([xmin2,xmax2,ymin2,ymax2])
hline(Pre1(:,:,j).limits(2),'b:',['VlimPre1= '  ,num2str(Pre1(:,:,j).limits(2))]);
vline(Pre1(:,:,j).limits(1),'b:',['PlimPre1= ' ,num2str(Pre1(:,:,j).limits(1))]);
hline(Pre2(:,:,j).limits(2),'g:',['VlimPre2= '  ,num2str(Pre2(:,:,j).limits(2))]);
vline(Pre2(:,:,j).limits(1),'g:',['PlimPre2= ' ,num2str(Pre2(:,:,j).limits(1))]);
vline(Post2(:,:,j).limits(1),'r:',['PlimPost= ',num2str(Post2(:,:,j).limits(1))]);

end
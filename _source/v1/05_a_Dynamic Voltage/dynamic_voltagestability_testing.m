%% VOLTAGE STABILITY
clear all; clc; close all

% 1. Chose between one of the two case studies to undestand the
% voltage stability function

%%%%%%%%%%%%%%%%% Artificial Simulation Resutls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load LOAD_ALL
%load LONNY_ALL
ta_pre=[1,5];
tb_pre=[6,10];
tc_pre=[11,15];

ta_post=[16,20];
tb_post=[21,25];
tc_post=[26,30];

Nb=2;
Nl=1;
Limits=[4,4];


%%%%%%%%%%%%%%%%% RTE Simulation Resutls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   load LONNY_ALL
% ta_pre=[30,200];
% tb_pre=[400,600];
% tc_pre=[800,900];
% 
% ta_post=[1000,1090];
% tb_post=[1100,1130];
% tc_post=[1140,1170];
% 
% %Pre_t=[30,200;400,600;800,900];
% %Pos_t=[1000,1090;1100,1130;1140,1170];
% Nb=10;
% %Nb=[5,10,22];   %Bus 10 Used for presentation for Oslo (Lonny)
% Nl=10; %Line 10 Used for presentation for Oslo (Lonny)
% Limits=[10,10];



%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%

% At least 3 set of simulation results has to be provided for 3 different loading
% levels (a=low, b=OK, c=limit), the signals required are: Time (t), Voltage (V), Active Power (P) and Reactive Power (Q) in cell array format
% where each element of the cell array, corresponds to a different simulation result for each loading level

t={t20,t60,t80}; % Time  
V={V20,V60,V80}; % Voltage at representative buses
P={P20,P60,P80}; % Active Power in representaive lines
Q={Q20,Q60,Q80}; % Reactive Power on representativr lines

t_pre = {ta_pre,tb_pre,tc_pre};   %Section of data used to estimate the PV curve in pre-contingency
t_post= {ta_post,tb_post,tc_post};%Section of data used to estimate the PV curve in post-contingency

Nb; %Number of bus or buses to analyze, i.e. Nb=[10] for single bus, Nb=[5,8,12] 3 buses analysis
Nl; %Number of line analized

Limits; % Limits=[Plim,Vlim] in (%),  Percent of Power from the Nose  

x01 = pi/8; % inital guess for angle deltas  
x02 = 1;    % intial guess for voltage  E
x03 = 0.03; % intial guess for reactance X

r01 = 1;  % Inital guess for alpha
r02 = 1;  % Initial guess for beta

x0=[x01;x02;x03];  
r0=[r01;r02];


%%

 [SBI,ABI,GBI,Pre,Post]=voltagestability(t,V,P,Q,t_pre,t_post,Nb,Nl,Limits,x0,r0);
 
 SBI
 ABI
 GBI
 
% Outputs
% SBI - Single Bus Index. Is a matrix and provides:  the distance in pre and post contingency 
%       for each loading level (low, ok and limit) to  the maximum loadability (Pmax) and voltage (Vlim)
%       limits for a selected bus(i).
% ABI - All Bus Index. Is a vector that provides the minimum distance among all buses for each 
%       loading level (low, ok and limit) to the loadability (Pmax) and the voltage (Vlim) limits.
%       The ABI index helps to identify if a loading level is violating a limit in pre or post contingency.
% GBI - Is a 2 element vector that provides the overall minimum distance to  the loadability (Pmax)
%       and the over all voltage (Vmax) limits respect to all buses. 
%       The GBI index indicates if a limit has been violated.
% Pre & Post - Structure containing the following information. 
%     .limits - row 1 Power limits, row 2 Voltage limits
%     .mean   - row 1, Power mean values for each loading level [Pa,Pb,Pb], where a=low, b=OK, c=limit levels 
%               row 2, Voltage mean values for each loading level [Va,Vb,Vc],where a=low, b=OK, c=limit levels 
%     .deltas - same information as SBI
%     .curve  - column 1 estimated Power, useful to reproduce nose curve
%               column 2 estimated Voltage, useful tor efproduce nose curve
%                    i.e. precontingency PV plot = plot(Pre.curve(:,1),Pre.curve(:,2))   
 
 
 
 return
%% Optional, PLotting Resulting Curves

qq=size(Nb,2); 

for j=1:qq,
pcent1=20/100;
pcent2=5/100;

Pm=Pre(:,:,j).limits(1)+(Pre(:,:,j).limits(1)*(Limits(1)/100));


xmin= Pre(:,:,j).mean(1,1)-( Pre(:,:,j).mean(1,1)*pcent1);
%xmax= Pre(:,:,j).limits(1)+(Pre(:,:,j).limits(1)*pcent2);
xmax= Pm+(Pm*pcent2);
ymin= Pre(:,:,j).limits(2)-(Pre(:,:,j).limits(2)*pcent1);
ymax= Pre(:,:,j).mean(2,1)+(Pre(:,:,j).mean(2,1)*pcent2);

figure; 
plot(Pre(:,:,j).mean(1,:),Pre(:,:,j).mean(2,:),'g*')
title('pre-contingency')
hold on;plot(Pre(:,:,j).curve(:,1),Pre(:,:,j).curve(:,2),'b-')
axis([xmin,xmax,ymin,ymax])

hline(Pre(:,:,j).limits(2),'b:',['Vlim= ',num2str(Pre(:,:,j).limits(2))]);
vline(Pre(:,:,j).limits(1),'b:',['Plim1= ',num2str(Pre(:,:,j).limits(1))]);
vline(Pre(:,:,j).mean(1,1),'k:',['Pa= ',num2str(Pre(:,:,j).mean(1,1))]);
vline(Pre(:,:,j).mean(1,2),'k:',['Pb= ',num2str(Pre(:,:,j).mean(1,2))]);
vline(Pre(:,:,j).mean(1,3),'k:',['Pc= ',num2str(Pre(:,:,j).mean(1,3))]);
hline(Pre(:,:,j).mean(2,1),'k:',['Va= ',num2str(Pre(:,:,j).mean(2,1))]);
hline(Pre(:,:,j).mean(2,2),'k:',['Vb= ',num2str(Pre(:,:,j).mean(2,2))]);
hline(Pre(:,:,j).mean(2,3),'k:',['Vc= ',num2str(Pre(:,:,j).mean(2,3))]);

 text(Pre(:,:,j).mean(1,1),Pre(:,:,j).mean(2,1),['\DeltaPa=',num2str(Pre(:,:,j).deltas(1,1))],'FontSize',10)
 text(Pre(:,:,j).mean(1,2),Pre(:,:,j).mean(2,1),['\DeltaPb=',num2str(Pre(:,:,j).deltas(1,2))],'FontSize',10)
 text(Pre(:,:,j).mean(1,3),Pre(:,:,j).mean(2,3),['\DeltaPc=',num2str(Pre(:,:,j).deltas(1,3))],'FontSize',10)
 text(Pre(:,:,j).mean(1,1),Pre(:,:,j).mean(2,2),['\DeltaVa=',num2str(Pre(:,:,j).deltas(2,1))],'FontSize',10)
 text(Pre(:,:,j).mean(1,2),Pre(:,:,j).mean(2,3),['\DeltaVb=',num2str(Pre(:,:,j).deltas(2,2))],'FontSize',10,'Color','k')
 text(Pre(:,:,j).mean(1,3),Pre(:,:,j).limits(2),['\DeltaVc=',num2str(Pre(:,:,j).deltas(2,3))],'FontSize',10)

figure; 
plot(Pre(:,:,j).mean(1,:),Pre(:,:,j).mean(2,:),'g*')
title('post-contingency')
hold on;plot(Pre(:,:,j).curve(:,1),Pre(:,:,j).curve(:,2),'b-')
hold on;plot(Post(:,:,j).curve(:,1),Post(:,:,j).curve(:,2),'r-')
hold on;plot(Post(:,:,j).mean(1,:),Post(:,:,j).mean(2,:),'c*')
axis([xmin,xmax,ymin,ymax])

hline(Pre(:,:,j).limits(2),'b:',['Vlim= ',num2str(Pre(:,:,j).limits(2))]);
vline(Pre(:,:,j).limits(1),'b:',['Plim1= ',num2str(Pre(:,:,j).limits(1))]);
vline(Post(:,:,j).limits(1),'r:',['Plim2= ',num2str(Post(:,:,j).limits(1))]);
vline(Pre(:,:,j).mean(1,1),'k:',['Pa= ',num2str(Post(:,:,j).mean(1,1))]);
vline(Pre(:,:,j).mean(1,2),'k:',['Pb= ',num2str(Post(:,:,j).mean(1,2))]);
vline(Pre(:,:,j).mean(1,3),'k:',['Pc= ',num2str(Post(:,:,j).mean(1,3))]);
hline(Post(:,:,j).mean(2,1),'r:',['Va= ',num2str(Post(:,:,j).mean(2,1))]);
hline(Post(:,:,j).mean(2,2),'r:',['Vb= ',num2str(Post(:,:,j).mean(2,2))]);
hline(Post(:,:,j).mean(2,3),'r:',['Vc= ',num2str(Post(:,:,j).mean(2,3))]);
 
 text(Post(:,:,j).mean(1,1),Post(:,:,j).mean(2,1),['\DeltaPa=',num2str(Post(:,:,j).deltas(1,1))],'FontSize',10,'Color','r')
 text(Post(:,:,j).mean(1,2),Post(:,:,j).mean(2,1),['\DeltaPb=',num2str(Post(:,:,j).deltas(1,2))],'FontSize',10,'Color','r')
 text(Post(:,:,j).mean(1,3),Post(:,:,j).mean(2,3),['\DeltaPc=',num2str(Post(:,:,j).deltas(1,3))],'FontSize',10,'Color','r')
 text(Post(:,:,j).mean(1,1),Post(:,:,j).mean(2,2),['\DeltaVa=',num2str(Post(:,:,j).deltas(2,1))],'FontSize',10,'Color','r')
 text(Post(:,:,j).mean(1,2),Post(:,:,j).mean(2,3),['\DeltaVb=',num2str(Post(:,:,j).deltas(2,2))],'FontSize',10,'Color','r')
 text(Post(:,:,j).mean(1,3),Post(:,:,j).limits(2),['\DeltaVc=',num2str(Post(:,:,j).deltas(2,3))],'FontSize',10,'Color','r') 
  
end

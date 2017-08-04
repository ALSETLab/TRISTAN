clear all;clc;close all

 %% Load time series
% 
 load SmallSignal
 
 

%% ERA & PRONY

Method =2;     % Method=1 (ERA analysis)
               % Method=2 (PRONY analysis recommended)
               
t;             % Time Vector
SigAn = LP1;   % Signals to analyze (Active power flows at critical lines)

t1=[1,35,25];  % Select window section t1(1)= Initial Time
                                     %  t1(2)= Final time
                                     %  t1(3)= From where to start ERA         
Pcent=100;     % Amount of signals to be considered, 100=100% for all signals,
f=[0.1,1];     % Range of frequencies to analyze damping of modes  f(1)=fmin, f(2)=fmax both in Hz
Nmodes=[];     % Number of Modes to estimate on each approach, Nmodes=[]; automatic selection, Nmodes=[2]; force Algorithm to choose 2 modes.
Damp=[0,5,10]; % Damping ratios, Calculate damping distance from each mode to the specified dampings, ie [0%, 5%, 10%]
flag1=1;       % if 1 show plots

[GMI,AMI,SMI,Poles,Freq,Ener]=smallsignal(Method,t,SigAn,t1,Pcent,f,Nmodes,Damp,flag1);

SMI
AMI
GMI

%Outputs
%
% Ener  Results of sorting signals "SigAn" according to energy
%       Ener.y  Sorted signals from t1(1) to t1(2) 
%       Ener.yh Sorted signals from t1(3) to t1(2)
%       Ener.t  Time vector from t1(1) to t1(2)
%       Ener.th Time vector from t1(2) to t1(2)
%
% Freq  Frequencies found from the FFT filter
%
% Poles Output of applying ERA or PRONY to set of signals "Ener.yh"
%
%      if ERA 
%      Poles.sys Estimated continous state-space model 
%
%      if PRONY 
%      Poles.sys Estimated "A" matrix with modes identified 
%
% SMI  Single Mode Indicator, Matrix. 
% AMI  All Mode Indicator, Vector.
% GMI  Global Mode Indicator, Gain.


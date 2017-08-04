clear all;clc;close all

%% Load time series

load Over_Load

%% POST-FAULT OVER LOADS

t;                 % Time vector
t1=[1,40,30];      % Section of Signal t1(1)= Initial Time
                   %                   t1(2)= Final time
                   %                   t1(3)= Number of samples prior t1(2)   

Signal = S1;       % Input signal, Aparent power (S) size t x N
nl     = size(S1,2); % Determine the size of Signal
wf_i   = ones(1,nl); % Uniform weights in all Lines 
p      = 3;            % Exponent used to magnify problems, value between 1 to "inf"
dev    = 10;           % Maximum Variation allowed from the nominal value (in %), i.e, 5 for 5%
flag   = 1;           %PLotting results


[f_x,f_scaled,f_noscaled]=static_overload(t,t1,Signal,wf_i,p,dev,flag);

f_x

% outputs
%      f_x        = Final index
%      f_scaled   = Indivial index of each Line with the exponent applyied
%      f_noscaled = Indivial index of each Line without exponent
    

%% eof
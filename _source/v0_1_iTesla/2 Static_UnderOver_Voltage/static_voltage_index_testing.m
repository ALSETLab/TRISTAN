% POST-FAULT UNDER/OVER VOLTAGE
% Sample script
clear all;clc;close all

%% Load time series
load OverUnder_Voltage

% POST-FAULT UNDER/OVER VOLTAGE
%% Input data and options into the index
t;                 % time vector

t1 =[1,900,50];    % Select section of Signal t1(1)= Initial Time
                   %                          t1(2)= Final time
                   %                          t1(3)= Number of samples before time t1(2)  
                   %                                 These make up the
                   %                                 window of data used to
                   %                                 compute the index.
                   
V    = BVm;        % Voltage mangitudes - Matrix of size t x N
nb   = size(BVm,2);% Determine the size of V
wv_i = ones(1,nb); % Uniform weights in all Buses 
                   % If information about each bus is provided, then a
                   % similar vector with the weights for each bus needs to
                   % be provided by the user
p    = 3;         % Exponent used to scale the index, value between 1 to "inf"
Dev  = 1;        % Maximum voltage variation allowed (in %) of the nominal value, i.e. 2 for 2%
Method = 0;        % Select the method for time series analysis
                   %     Method=0 average value from analyzed time window (recomendend)
                   %     Method=1 largest value from analyzed time window
flag =1;                      
% Compute the index
[v_x,v_scaled,v_noscaled] = static_voltage(t,t1,BVm,wv_i,p,Dev,Method,flag);
% Display the index in the Command Window
v_x

% outputs
%      v_x        = Final index
%      v_scaled   = Indivial index of each bus with the exponent applyied
%      v_noscaled = Indivial index of each bus without exponent
%% eof
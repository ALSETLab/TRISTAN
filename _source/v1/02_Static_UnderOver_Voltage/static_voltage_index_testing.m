% POST-FAULT UNDER/OVER VOLTAGE

% Sample script
clear all;clc;close all

%% Load time series
load OverUnder_Voltage

% POST-FAULT UNDER/OVER VOLTAGE
%% Input data and options into the index
t;               % time vector
V    = BVm;      % Voltage magnitudes - Matrix of size t x N
p    = 3;        % Exponent used to scale the index, value between 1 to "inf"
d    = 1;        % Maximum voltage variation allowed (in %) of the nominal value, i.e. 2 for 2%
simTime=max(t);  % Total simulation time of the given data
deltime=25;      %Time period of each interval


% Compute the index
[v_x,i2h,F,f,G,g,Steady_State_buses] = static_voltage(t,V,p,d,deltime,simTime);

% Display the index in the Command Window
v_x
%%
plot(t,V);
if SS_buses~=0
figure;plot(t,V(:,SS_buses(:,:)))
end
% outputs
%  v_x            = Final index
%  F              = voltage level of each bus as per the index calculation
%  f              = Bus numbers which violated voltage limits
%  Steady_State_buses = Bus numbers whose voltages reached steady state
%  G              = Rate of change of voltage variation of all buses
%  g              = Buses in high voltage variation
%% eof
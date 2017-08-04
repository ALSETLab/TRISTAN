clear all;clc;close all

%% Load time series

load Over_Load
%load Overload-1
%load Overload-2
%load Overload-3
%% POST-FAULT OVER LOADS

t;             % Time vector
S  = S1;       % Input signal, Aparent power (S) size t x N
p  = 3;        % Exponent used to magnify problems, value between 1 to "inf"
d  = 10;       % Maximum Variation allowed from the nominal value (in %), e.g, 10 for 10%
simTime=max(t);% Total simulation time of the given data
Faultime=10;   % Fault occuring time in the data
deltime=10;    %Time period of each sampled interval

figure;plot(t,S);

[f_x,i2h,ridh,F,f,G,g,steady_state_lines] = static_overload(t,S,p,d,deltime,simTime);

f_x
% outputs
%      f_x         = Final index
%      F           = Level of loading on each line
%      f           = Overloaded lines
%      G           = State of the lines
%      g           = The lines in which the power flow variation are high
%steady_state_lines= The lines which reached steady state
%% eof
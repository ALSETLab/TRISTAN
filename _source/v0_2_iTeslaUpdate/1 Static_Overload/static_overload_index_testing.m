clear all;clc;close all

%% Load time series

load Over_Load

%% POST-FAULT OVER LOADS

t;                 % Time vector
S  = S1;       % Input signal, Aparent power (S) size t x N
p  = 3;        % Exponent used to magnify problems, value between 1 to "inf"
d  = 10;       % Maximum Variation allowed from the nominal value (in %), e.g, 10 for 10%


 f_x = static_overload(t,S,p,d);

f_x

% outputs
%      f_x        = Final index
  

%% eof
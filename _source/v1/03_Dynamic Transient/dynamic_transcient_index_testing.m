clear all;clc;close all

%% Load time series

load Transient

%% TRANSIENT STABILITY

t;                 % Time vector size t x 1                
delta;  % Angle of the machines size t x N;
M;      % Two time interia (H) of the machines (M=2*H), size 1 x N

[J]=dynamic_transient(t,delta,M)


% outputs
%      J       = Final index



clear all;clc;close all

%% Load time series

load Transient

%% TRANSIENT STABILITY

t;                 % Time vector size t x 1
t1=[1,63];     % Section of Signal t1(1)= Initial Time in seconds
                   %                   t1(2)= Final time in seconds
                    
delta;  % Angle of the machines size t x N;
M;      % Two time interia (H) of the machines (M=2*H), size 1 x N

[J]=dynamic_transient(t,t1,delta,M)


% outputs
%      J       = Final index



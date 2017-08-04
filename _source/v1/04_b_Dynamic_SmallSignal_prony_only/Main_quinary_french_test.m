clear all;clc;close all
addpath 'H:\KTH-Postdoc 2013\iTesla\WP4\Indexes_review\20130115_1845_SN2_FR0'

%% Load time series
 
 load sim_fault_0_smallsignal.mat
% load sim_fault_1_smallsignal.mat
% load sim_fault_2_smallsignal.mat
% load sim_fault_3_smallsignal.mat
% load sim_fault_4_smallsignal.mat
% load sim_fault_5_smallsignal.mat
% load sim_fault_6_smallsignal.mat
% load sim_fault_7_smallsignal.mat
 
figure;plot(t,LP1);axis tight 
return
 %% Set input parameters and call main function

time   = t;      % Time vector of order Nx1 with variable time step
signal  = LP1; % Active power flow on relevant lines of order NxM, where M is the number of signals

%step_min = 0.04;  % Minimum step size (full details in future documentation)
%var_min = 0.1%0.3;    % [0-1] Filter signals with varianza smaller than "var_min" 
%f = [0.1,2.5];    % [fmin,fmax] Range of  frequencies of interest in Hz
%d =[0,5,10];      % Damping values in percent (%) where the index distances will be calculated
%Nm  =10;          % Number of modes used in Prony to reconstruct the signals

[ss,y0,G,det]=sssi(signal,time,step_min,var_min,f,d,Nm);

ss.smi
ss.ami
ss.gmi

%% Outputs
%
% ss is a structure with the actual small-signal-stability index
%       ss.smi single mode indicator
%       ss.ami all modes indicator
%       ss.gmi global mode indicator
%
% y0 is a structure with the section of signals that where used as input in Prony  
%       y0.t   time vector 
%       y0.y   signals
%
% G is a structure with Prony output 
%       G.Poles detected modes 
%       G.Res   residue of each mode
%       G.K     gain of each mode
%       G.that  time vector
%       G.yhat  re-constructed signal
%
% det is a structure summarizing information
%       det.l  Complex value of detected modes in the range of frequency
%       det.f  Frequency in Hz of modes
%       det.d Damping in percent of modes
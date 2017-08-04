function [lamda,model]=pronyiTesla(t,y,XX,tstart,tend,T,shift,tstplot,tedplot,plotFlag);
%[model]=pronyanalysis(t,y,n,tstart,tend,T,shift,tstplot,tedplot,plotFlag);
%Prony analysis program for fitting to a ringdown
% Inputs:
%   t = time vector (does not need to be equally spaced) -- column vector
%       t(1) is assumed to be 0; if not, all time variables are shifted by t(1).
%   y = ringdown matrix of order N by NChan corresponding to t, each column is a different signal
%   n = order of the estimated model
%   tstart = row vector of order 1 by NChan.  tstart(k) is the starting
%            time for analysis for y(:,k).
%   tend = ending times for analysis (same dimension as tstart)
%   shift = flag; if shift = 1, residues are shifted to t = 0.
%           If the data is noisy, reccomend shift = 0.
%   T = sample period for analysis
%   tsttplot = starting time for model simulation
%   tedplot = ending time for model simulation
%   plotFlag = if = 1, plot results
% Outputs (structured array):
%   model.Poles = ringdown pole matrix -- column k is for column k of y
%   model.Res = ringdown residue matrix
%   model.K = ringdown offset row vector
%   model.yhat = model signal matrix
%   model.that = time vector for yhat (starts at tstart)
%
% NOTE: It is reccomended that the N/(1+NChan) < 200, where
%       N = total number of data points in y to be analyzed, and
%       NChan = number of columns in y. 

% Written by D. Trudnowski, 1999.
% Last update, D. Trudnowski, 2005.

% 1.0 Setup

%Basic error checks

if (size(y,2)~=size(tstart,2))|(size(y,2)~=size(tend,2)); error('Dimension error in y, tstart, tend'); end;
if (size(tstart,1)~=1)|(size(tend,1)~=1); error('Dimension error in tstart or tend'); end;
if (T<=0)|(tstplot>=tedplot)|(max(t)<max(tend'))|(min(t)>min(tstart')); error('data error'); end;
if (size(t,1) ~= size(y,1))|(size(t,2)~=1); error('Dimension error in y or t'); end;

%Data parameters
NChannels = size(y,2);

%Shift time parameters
tstart = tstart-t(1);
tend = tend-t(1);
tstplot = tstplot-t(1);
tedplot = tedplot-t(1);
t = t-t(1);

%Set up analysis data and calculate offset
tanalysis = T*[0:1:ceil(max(tend')/T)+1]';
Nstart = zeros(1,NChannels);
Nend = zeros(1,NChannels);
yanal = zeros(length(tanalysis),NChannels);
K = zeros(1,NChannels);
for k=1:NChannels
    Nstart(k) = floor(tstart(k)/T)+1;
    Nend(k) = ceil(tend(k)/T)+1;
    yanal(:,k) = spline(t,y(:,k),tanalysis);
    K(1,k) = mean(yanal(Nstart(k):Nend(k),k));
    yanal(:,k) = yanal(:,k)-K(k);
end
clear k

%Find model order
NdataPoints = Nend-Nstart+1; %Number of data points used for analysis on each channel
Ntotal = sum(NdataPoints'); %Total number of data points used for Prony analysis
%nOrder = round(Ntotal/(1+NChannels))-10; %Order of Linear Prediction
nOrder=XX;
%if (nOrder>50); nOrder=50; end; %Limit order to avoid numerical problems

% 2.0 Solve Linear Prediction

%Build matrix and vector
for k=1:NChannels
    Ym = zeros(NdataPoints(k)-nOrder,nOrder);
    for kk=1:nOrder
        Ym(:,kk) = yanal(Nstart(k)+nOrder-kk:Nstart(k)-kk+NdataPoints(k)-1,k);
    end
    yv = yanal(Nstart(k)+nOrder:Nstart(k)+NdataPoints(k)-1,k);
    if k==1;
        Ymatrix = Ym;
        yvector = yv;
    else
        Ymatrix = [Ymatrix;Ym]; %Cancatinate the channels
        yvector = [yvector;yv];
    end
end
clear Ym yv k kk
% acoef = pinv(Ymatrix)*yvector; %characteristic eqn. Pseudo inverse using SVD
acoef = Ymatrix\yvector; %This provides a much more accurate solution than the above line

%Find poles
zPoles = roots([1;-acoef]);
sPoles = log(zPoles)/T;



% 3.0 Solve for residues
Res = zeros(nOrder,NChannels);
for k=1:NChannels
    ZMatrix = zeros(NdataPoints(k),nOrder);
    for kk=1:NdataPoints(k)
        ZMatrix(kk,:) = (zPoles.').^(kk-1);
    end
    Res(:,k) = ZMatrix\yanal(Nstart(k):Nend(k),k);
    if shift==1; 
        Res(:,k) = Res(:,k).*(zPoles.^(-Nstart(k)+1)); %Shift residues to time 0
    end
end
clear k kk

% 4.0 Re-order using energy
P = zeros(nOrder,NChannels);
R = zeros(size(Res));
for k=1:NChannels
    clear E
    for kk=1:nOrder
        if abs(real(sPoles(kk)))<1e-8
            E(kk)=Res(kk,k)^2*(tend(k)-tstart(k));
        else
            E(kk)=(Res(kk,k)^2/(2*real(sPoles(kk))))*(exp(2*real(sPoles(kk))*(tend(k)-tstart(k)))-1);
        end
    end
    E=E(:);
    [x,ii]=sort(E);
    R(:,k)=Res(ii,k); 
    P(:,k)=sPoles(ii);
    M=[length(ii):-1:1]';
    R(:,k)=R(M,k); 
    P(:,k)=P(M,k);
end
clear M k x ii E

%Simulate;
that = [tstplot:T:tedplot]';
yhat = zeros(length(that),NChannels);
for k=1:NChannels
    yhat(:,k) = K(k).*ones(size(that));
    for kk=1:length(that);
        for n=1:nOrder;
            yhat(kk,k) = yhat(kk,k) + R(n,k)*exp(P(n,k)*that(kk));
        end
    end;
end;
yhat=real(yhat);

% 5.0 Place output in structured array
model.Poles = P;
model.Res = R;
model.K = K;
model.that = that;
model.yhat = yhat;

lamda=sPoles;



% 6.0 Plot results
if plotFlag==1
       figure
       hold on
        h1 = plot(t,y);
        h2 = plot(that,yhat,'--*');axis tight
        hold off
        xlabel('Time (sec)')
        %legend('Actual','Prony')
end

%Original plots
% if plotFlag==1
%     for k=1:size(yhat,2)
%         figure
%         n = find(t(t>=tstplot&t<=tedplot));
%         hold on
%         h1 = plot(t(n),y(n,k),'k','linewidth',1.5);
%         h2 = plot(that,yhat(:,k),'r--','linewidth',1.5);
%         hold off
%         xlabel('Time (sec)')
%         legend('Actual','Prony')
%     end
% end
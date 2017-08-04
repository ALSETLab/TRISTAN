function   [ymm,tss,ymm_sec,tsec]=energy_sort(t,t1,BVm,Pcent,flag1)
%
% [ymm,ymm_sec,tsec]=energy_sort(t,t1,Signals,Energy%,flag1)
%
%
%
%

%% Pre-processing the signals: Define the window of time to be analyzed

%  t1 =[tnom, tend, dt]   with the following information:
%
%   tnom -> Initial time to be analized (sec) 
%   tend -> End time to be analized (sec)
%   dt   ->Time to start analyzing small-signal  
  

% The following section extracts from the time vector "t",
% the approximated values of "tnom" and  "tend" defined in t1

tnom0 = find(t>=t1(1)*0.95 & t<=t1(1)*1.05 );
    if size(tnom0,1)>=3
            tnom0 = find(t>=t1(1)*0.999 & t<=t1(1)*1.001 );
    end
    if size(tnom0,1)==0
            tnom0 = find(t>=t1(1)*0.9 & t<=t1(1)*1.1 );
    end
   
tend0 = find(t>=t1(2)*0.95 & t<=t1(2)*1.05);
    if size(tend0,1)>=3
            tend0 = find(t>=t1(2)*0.999 & t<=t1(2)*1.001 );
    end
    if size(tend0,1)==0
            tend0 = find(t>=t1(2)*0.9 & t<=t1(2)*1.1 );
    end

tnom = tnom0(1);
tend = tend0(end);
tdel = t(tnom:tend);
  %dt = t1(3);

  
  
%% Selecting voltage magnitudes with largest oscillation content 
nb       = size(BVm,2);
tss      = t(tnom:1:tend); % Time window to analyze
tdss     = [tnom:1:tend];  % 
y_invest = BVm(tdss,:);    % Voltage magnitudes using selected interval


tosc0 = find(tss>=t1(3)*0.95 & tss<=t1(3)*1.05);
    if size(tosc0,1)>=3
            tosc0 = find(tss>=t1(3)*0.999 & tss<=t1(3)*1.001 );
    end
    if size(tosc0,1)==0
            tosc0 = find(tss>=t1(3)*0.9 & tss<=t1(3)*1.1 );
    end   
  dt=tosc0(end);

% Calculating the Energy of the signals
for i=1:nb
     yme(:,i)=y_invest(:,i)-mean(y_invest(:,i)); % Detrended signals
     %ydtr(:,i)=dtrend(y_invest(:,i)); % Detrended using matlab function
     En_sigs(i)=sqrt(sum(abs(yme(:,i).^2)));     % Energy  
end


[Ener,Ind] = sort(En_sigs','descend'); % Sort energy from largest to smallest

Pcent=100-Pcent;
Emax=max(Ener);
Emin=(max(Ener)*Pcent)/100;

kk=1;
for i=1:nb
if Ener(i)>= Emin
    EnerM(kk,1)=Ener(i);
    IndM(kk,1)=Ind(i);
    kk=kk+1;
end
end

 y_dtrend=yme(:,IndM);
 ymm=y_dtrend;  
%  ymm_sec=ymm((end-dt:end),:);
%  tsec=tss(end-dt:end);
%ymm_sec=ymm((dt:end),:);
ymm_sec=dtrend(ymm((dt:end),:));
tsec=tss(dt:end);




% y_dtrend=yme(:,Ind(1:end));
% ymm=y_dtrend;  
% ymm_sec=ymm((end-dt:end),:);
% tsec=tss(end-dt:end);

%   keyboard

%%

if flag1==1
%figure;plot(tss,BVm(tdss,Ind(1:end)));axis tight % Plot the most significant signals
figure;plot(tss,ymm);axis tight % Plot the most significant signals detrended
title(['Signals sorted according to Enery, more than ',num2str(Pcent),' %'])
%figure;plot(tss(end-dt:end,1),ymm(end-dt:end,:));axis tight
figure;plot(tss(dt:end,1),ymm_sec);axis tight
title(['Selected Window to analyze: ',num2str(dt),' samplings before ', num2str(t(tend)),' sec'])
end
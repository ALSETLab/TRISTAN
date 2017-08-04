
function [Over_S,i2h,ridh]=static_overload(t,t1,S1,wf_i,p_1,dev0, flag1);
%% Post-fault Over Load Index
% This index is used to observer if the post-fault flows surpass the network capacity,
% by monitoring the power flows through the transmission lines right after an outage occurs.
%
%
% [Sindex,SingleIndex]=static_overload(t,t1,Signal,w,p,dev);
%
%  t  - Time vector
%  t1 - [tnom, tend, dt] 
%        tnom - Initial time to be analized (sec) 
%        tend - Final time to be analized (sec)
%          dt - Number of sampling times before the final time, to be analyzed
%  Signal - Signal to analyze, Apparent power flow 
%  w - Weights for all buses, vector of ones for unitary weight
%  p - Exponent
%  dev - Deviation allowed of power flow from nominal value, i.e 10 for +10%    
% 
%                                   Version 1.1

dev0=dev0/100;
%% Pre-processing the signals: Define the window of time to be analyzed

%  t1 =[tnom, tend, dt]   with the following information:
%
%   tnom -> Initial time to be analized (sec) 
%   tend -> End time to be analized (sec)
%   dt   -> Number of samplings to be analyzed, i.e. dt=3, t_analyzed = [t(tend-3),t(tend-2),t(tend-1),t(tend)]  
  

% The following section extracts from the time vector "t",
% the approximated values of "tnom" and  "tend" defined in t1

tnom0 = find(t>=t1(1)*0.95 & t<=t1(1)*1.05 );
    if size(tnom0,1)==0
            tnom0 = find(t>=t1(1)*0.9 & t<=t1(1)*1.1 );
   end
tend0 = find(t>=t1(2)*0.95 & t<=t1(2)*1.05);
    if size(tend0,1)==0
            tend0 = find(t>=t1(2)*0.9 & t<=t1(2)*1.1 );
   end


tnom = tnom0(1);
tend = tend0(end);
tdel = t(tnom:tend);
  dt = t1(3);



%% OVERLOAD 

nl=size(S1,2);
Snom = S1(tnom,:);        %  nominal value (pre-contingency)
Smax = Snom*(1+dev0);     %  Over power flow limits 
Spost = S1(tend-dt:tend,:); 
zer = find(Smax<=0.001);

if zer>1
     Smax(zer)=0.1;
end


Smean = mean(Spost);

k=1; over_line=[];
for i=1:nl
indxs_loc=(wf_i(i)*((abs(Smean(1,i))/abs(Smax(1,i)))^p_1));
index_red(i)=indxs_loc;


if indxs_loc>=1
   indx(i)=indxs_loc;
   over_line(k)=i; % indices of lines violating limits
   k=k+1;
else
   indx(i)=1;
end


rid(i,1) = index_red(i).^(1/p_1);
if rid(i)<=1
    ridh(i,1)=1;
    i2h(i,1)=1;
    else               
   ridh(i,1)= rid(i); 
   i2h(i,1)=index_red(i);
end

end


Over_S=sum(i2h)/nl;

%%
if flag1==1

if Over_S>=1
    nvl=size(over_line,2); %Number of lines violating the limits
    Sextra_pu= indx(over_line)/max(indx(over_line));
    for i=1:nvl
        if  indx(over_line(i)).^(1/p_1)<=2
    indx_over_per(i,1)=((indx(over_line(i)).^(1/p_1))-1)*100;
        else
             indx_over_per(i,1)=((indx(over_line(i)).^(1/p_1)))*100;
        end
    end
    
 if nvl>=1
    
figure;bar(Smean(1,over_line),'g');hold on;bar(Smax(1,over_line),'y'); legend('S_{mean}','S_{max_{limit}}')
xlabel('Line Number');ylabel('MVA'); title(['Lines violating the limits=',num2str(nvl)]);
set(gca,'XTick',[1:nvl])
set(gca,'XTickLabel',over_line)


figure;bar(Sextra_pu,'r')
rr=find(indx==max(indx));
title(['Line ',num2str(rr),' with maximum violation of ',num2str(max(indx_over_per)),'% more'])
set(gca,'XTick',[1:nvl])
set(gca,'XTickLabel',over_line)


%figure;plot(t,S1(:,rr));axis tight;title(['Line number ', num2str(rr),' with maximum violation'])
figure;plot(tdel,S1(tnom:tend,rr));axis tight;title(['Line number ', num2str(rr),' with maximum violation'])
hold on; plot(t(tend-dt:tend,1),S1(tend-dt:tend,rr),'linewidth',3)
%hold on; hline(Smax(:,rr),'k--');hold on; hline(mean(Spost(:,rr)),'r')
hold on; horline(tdel,Smax(:,rr),'k:');hold on; horline(tdel,mean(Spost(:,rr)),'r')
xlabel('Time (sec)'); ylabel('MVA')
 end

end

end


 
end
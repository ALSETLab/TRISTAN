
function [Over_S,i2h,ridh,F,f,G,g,steady_state_lines]=static_overload(t,S1,p_1,dev0,deltime,simTime)
%% Post-fault Over Load Index
% This index is used to observe if the post-fault flows surpass the network capacity,
% by monitoring the power flows through the transmission lines right after an outage occurs.
%
%
% [Sindex]=static_overload(t,t1,Signal,w,p,d);
%
% INPUTS
%  t - Time vector
%  S - Apparent power flow 
%  p - Exponent
%  d - Deviation allowed of power flow from nominal value, e.g. 10 for +10%    
% 
%  OUTPUTS
%  Sindex -Overload index 
%
%                                   Version 1.2


%% OVERLOAD 
pre0  = 5;
post0 = 100;

dev0  = dev0/100;  
nl    = size(S1,2);  
wf_i  = ones(1,nl);         % Uniform weights in all Lines is assumed 
Snom  = mean(S1(1:pre0,:)); % Nominal value (pre-contingency)
Smax  = Snom*(1+dev0);      % Over power flow limits 
Spost = S1(end-post0:end,:); 
Smean = mean(Spost);

% Remove lines out of service
xxx=[];k=1;
for i=1:nl
 if abs(Smean(1,i))<1e-2
     xxx(k,1)=i;
     S1(:,i)=0;
     Smean(1,i)=0;
     Snom(1,i)=0;
     k=k+1;
 end    
end
clear k

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
for i=1:size(rid)
Fx(i)=rid(i);
FFx(i)=i;
end
F=[FFx;Fx];
F=F';
Overloaded_lines=over_line';
for i=1:size(Overloaded_lines)
fx(i)=ridh(Overloaded_lines(i));
ffx(i)=Overloaded_lines(i);
end
f=[ffx;fx];
f=f';
Over_S=sum(i2h)/nl;

%% Plotting results (if required)
% if Over_S>=1
%     nvl=size(over_line,2); %Number of lines violating the limits
%     Sextra_pu= indx(over_line)/max(indx(over_line));
%     for i=1:nvl
%         if  indx(over_line(i)).^(1/p_1)<=2
%     indx_over_per(i,1)=((indx(over_line(i)).^(1/p_1))-1)*100;
%         else
%              indx_over_per(i,1)=((indx(over_line(i)).^(1/p_1)))*100;
%         end
%     end
%     
%  if nvl>=1
%     
% figure;bar(Smean(1,over_line),'g');hold on;bar(Smax(1,over_line),'y'); legend('S_{mean}','S_{max_{limit}}')
% xlabel('Line Number');ylabel('MVA'); title(['Lines violating the limits=',num2str(nvl)]);
% set(gca,'XTick',[1:nvl])
% set(gca,'XTickLabel',over_line)
% 
% figure;bar(Sextra_pu,'r')
% rr=find(indx==max(indx));
% title(['Line ',num2str(rr),' with maximum violation of ',num2str(max(indx_over_per)),'% more'])
% set(gca,'XTick',[1:nvl])
% set(gca,'XTickLabel',over_line)
% 
% figure;plot(t,S1(:,rr));axis tight;title(['Line number ', num2str(rr),' with maximum violation'])
% hold on; plot(t(end-post0:end),S1(end-post0:end,rr),'linewidth',3)
% hold on; horline(t,Smax(:,rr),'k:');hold on; horline(t,mean(Spost(:,rr)),'r')
% xlabel('Time (sec)'); ylabel('MVA')
% end
% 
% end
%%
t1=simTime-3*deltime;
t2=simTime-2*deltime;
t3=simTime-deltime;
t4=simTime;
S=S1;
lines=size(S);   

% Sampling three equal intervals towards the end of the simulation for 'k-th' line
for k=1:lines(2)
count1=0;count2=0;count3=0;
for i=1:length(t)    
if(t(i)>t1 && t(i)<=t2)
count1=count1+1;
SS1(count1,k)=S(i,k);
tt1(count1,k)=t(i);
else if (t(i)>t2 && t(i)<=t3)
count2=count2+1;
SS2(count2,k)=S(i,k);
tt2(count2,k)=t(i);
else if (t(i)>t3 && t(i)<=t4)
count3=count3+1;
SS3(count3,k)=S(i,k);
tt3(count3,k)=t(i);
    end
    end
end
end

% Calculating the slopes of the intervals for 'k-th' line
for i=2:count1
Slope1(i-1,k)=(SS1(i,k)-SS1(i-1,k))/(tt1(i,k)-tt1(i-1,k));
if Slope1(i-1,k)== Inf || Slope1(i-1,k)==-Inf 
Slope1(i-1,k)=0;
end
end

for i=2:count2
Slope2(i-1,k)=(SS2(i,k)-SS2(i-1,k))/(tt2(i,k)-tt2(i-1,k));
if Slope2(i-1,k)== Inf || Slope2(i-1,k)==-Inf
Slope2(i-1,k)=0;
end
end

for i=2:count3
Slope3(i-1,k)=(SS3(i,k)-SS3(i-1,k))/(tt3(i,k)-tt3(i-1,k));
end

% Calculation of mean slope of the interval for 'k-th' line
 mean_slope(1,k)=mean(Slope1(:,k));
 mean_slope(2,k)=mean(Slope2(:,k));
 mean_slope(3,k)=mean(Slope3(:,k));
 
% Calculation of slope variation between the intervals for 'k-th' line 
 slope_change1(k,1)=(mean_slope(2,k)-mean_slope(1,k));
 slope_change1(k,2)=(mean_slope(3,k)-mean_slope(2,k));

% Calculation of the difference between the slope variation for 'k-th' line  
variation(k,1)=(abs(slope_change1(k,2)-slope_change1(k,1)));
end

% Assigning all the lines to an array
 for m=1:lines(2)
 mm(m,1)=m;
 end

G=[mm variation]; % Matrix with slope variation of all lines

c1=0;   %counting variable
% loop to find the lines in which the power change/variation is high
for i=1:lines(2)
if (variation(i)>2.5)
    c1=c1+1;
    GG(c1,1)=variation(i);
    gg(c1,1)=mm(i);
end
end
% if all lines reaches steady state then setting g to zero
if c1>0
g=[gg,GG];
else
g=[0,0];
end
[mm slope_change1];
[mm variation];
% figure;plot(slope_change1(:,1),'-ob');
% hold on;plot(slope_change1(:,2),'-or');
% Finding the lines in which reached steady state
Slines=0;   % variable to count the lines which reached steady state
for i=1:lines(2)
if (slope_change1(i,1)>0 && slope_change1(i,2)<0 || abs(slope_change1(i,2))<=abs(slope_change1(i,1)) && variation(i,1)<0.01)
    Slines=Slines+1;
    steady_state_lines(Slines,1)=i;
end
end
if Slines==0
   steady_state_lines=0;
end

end
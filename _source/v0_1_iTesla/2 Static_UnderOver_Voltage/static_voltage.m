function [V_index,i2h,ridh]=static_voltage(t,t1,BVm,wv_i,p_2,dev,flag1,flag2);
%% Post-fault Under/Over Voltage Index 
% Index for voltage of the transmission network right after an outage has occurred,
% the index indicates if the voltage (BVm), surpass the minimum or
% maximum limits of the operational standards
%
%
% [Vindex,single_indexscaled,single_index]=satatic_voltage(t,t1,V,w,p,dev,Method);
%
%  t  - Time vector
%  t1 - [tnom, tend, dt] 
%        tnom - Initial time to be analized (sec) 
%        tend - Final time to be analized (sec)
%          dt - Number of sampling times before the final time, to be analyzed
%  V - Voltages vector 
%  w - Weights for all buses, vector of ones for unitary weight
%  p - Exponent
%  dev - Deviation allowed of voltage from nominal value, i.e 2 for +-2%    
%  Method - 0 to use the average value from analized window
%          1 to take the largest value from analyzed window
%
%                       Version 1.0
% 
dev=dev/100;

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

%%  Remove any buses that are out of service

Vfault = [];
Vfault = find(BVm(tend,:)==0); %Line out of service
nbu=0;
if size(Vfault,2)>0
    nbu=size(Vfault,2);
    BVm=BVm(:,[1:Vfault-1,Vfault+1:end]);   
%     figure(101);plot(t,BVm);title(['Voltage magnitudes, without the faulty line ',num2str(Vfault)])
%     title(['Voltage magnitudes, fault in line ',num2str(Vfault)])
end

%% The following parameters are used to calculate the index:

tpost = t(tend-dt:tend);      % Period of time analyzed
Vnom =  BVm(tnom,:);          % Nominal  value (pre-contingency)
Vpost = BVm(tend-dt:tend,:);  % Voltage values post-contingency to be analyzed
Vmin =  Vnom*(1-dev);         % Under voltage limits 
Vmax =  Vnom*(1+dev);         % Over voltage limits 
  nb =  size(Vnom,2);         % Number of buses 


delV   = (Vmax - Vmin)/2;     %Equation 
Vmean  = mean(Vpost);         %Equation


%% Calculation of the actual index

for i=1:nb-nbu

    
if flag1==1;
    
    % Index equaiton based on several points
    indxs =wv_i(1,i)*((abs(Vnom(i)-Vpost(:,i))./delV(i)).^p_2)';
else
    
    % Index equation base on average value    
    indxs =wv_i(1,i)*((abs(Vnom(i)-Vmean(:,i))./delV(i)).^p_2)';
end



% Applying considerations to understand the value of the index
i2(i,1)=max(indxs);


rid(i,1) = i2(i,1).^(1/p_2);

if rid(i)<=1
    ridh(i,1)=1;
    i2h(i,1)=1;
else if (rid(i)>=1) & (rid(i)<=5*(dev*100))
     ridh(i,1)=rid(i);
     i2h(i,1)=i2(i);
    else               
   ridh(i,1)= rid(i);
   i2h(i,1)=i2(i);
    end
end


end


i2_tot=sum(i2);

%%%%%%%%%%%% Under/Over Voltage Index %%%%%%%%%
V_index = sum(i2h)/(nb-nbu);  % The actual index normalized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Ploting results
 r1=find(i2>=1);   % index of lines violating limits
 nvl= size(r1,1);  % number of lines violating limits

 r2=(i2).^(1/p_2);
 %nvlsevere=size(r20,1); % number of lines violating the limits more than "n" times the variation (delV) allowed.
nvlsevere=0; 



if flag2 ==1

 if nvl>0
figure;plot(t(tend-dt:tend,1),BVm(tend-dt:tend,r1));axis tight
title(['Lines violating limits = ', num2str(nvl),'   Index=',num2str(V_index)])
if nvlsevere>0
%figure;plot(t,BVm(:,r20));axis tight
figure;plot(tdel,BVm(tnom:tend,r20));axis tight
title(['Lines (',num2str(nvlsevere), ') violating limits ',num2str(dev*100),' %'])
end
a= find(i2(r1)==max(i2(r1)));
%a= a(1);
mm=r1(1);
%figure;plot(t,BVm(:,r1(a)));axis tight

if flag1==1;

figure;plot(tdel,BVm(tnom:tend,r1(mm)));axis tight
hold on; plot(t(tend-dt:tend,1),BVm(tend-dt:tend,r1(mm)),'linewidth',3)
%hold on; hline(Vmin(r1(mm)),'k--');hline(Vmax(r1(mm)),'k--')
hold on; horline(tdel,Vmin(r1(mm)),'k:');horline(tdel,Vmax(r1(mm)),'k:')
title(['Line with max violation ', num2str(r1(mm)),' from the ',num2str(nvl),' lines violating the limits (',num2str(dev*100),' %)'])

else
  
figure;plot(tdel,BVm(tnom:tend,mm));axis tight
hold on; plot(t(tend-dt:tend,1),BVm(tend-dt:tend,mm),'linewidth',3)
%hold on; hline(Vmin(mm),'k--');hline(Vmax(mm),'k--')
%hold on; hline(mean(Vpost(:,mm)),'r')
hold on; horline(tdel,Vmin(mm),'k:');horline(tdel,Vmax(mm),'k:')
hold on; horline(tdel,mean(Vpost(:,mm)),'r')
title(['Line ', num2str(mm),' with maximum violation,  v=',num2str(V_index)])
xlabel('Time (sec)')
ylabel('KV')

end

 end
 
 end


 
 end


 
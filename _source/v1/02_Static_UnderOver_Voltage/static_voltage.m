function [V_index,i2h,ridh,Vol_violation,G,g,steady_state_bus]=static_voltage(t,BVm,p_2,dev,deltime,simTime)
%% Post-fault Under/Over Voltage Index 
% Index for voltage of the transmission network right after an outage has occurred,
% the index indicates if the voltage (V), surpass the limits of the operational standards
%
%
% [Vindex]=satatic_voltage(t,V,p,d);
%
% INPUTS
%  t - Time vector
%  V - Voltages vector 
%  p - Exponent
%  d - Deviation allowed of voltage from nominal value, e.g. 2 for +-2%    
%
% OUTPUTS
%  Vindex -Under/Over voltage index 
%
%
%                       Version 1.1
% 

dev=dev/100;


%% The following parameters are used to calculate the index:
 pre0  = 5;
 post0 = 100;


tpost = t(end-post0:end);     % Period of time analyzed
Vnom =  mean(BVm(1:pre0,:));  % Nominal  value (pre-contingency)
Vpost = BVm(end-post0:end,:); % Voltage values post-contingency to be analyzed
Vmin = Vnom*(1-dev);         % Under voltage limits 
Vmax = Vnom*(1+dev);         % Over voltage limits 
nb =  size(Vnom,2);           % Number of buses 
wv_i = ones(1,nb);            % Uniform weights in all Buses  


delV   = (Vmax - Vmin)/2;     %Equation 
Vmean  = mean(Vpost);         %Equation


% Remove bus out of service
xxx=[];k=1;
for i=1:nb
 if abs(Vmean(1,i))<1e-2
     xxx(k,1)=i;
     BVm(:,i)=0;
     Vmean(1,i)=0;
     Vnom(1,i)=0;
     Vmin(1,i)=0;
     Vmax(1,i)=0;
     wv_i(1,i)=0;
     k=k+1;
 end    
end
clear k

over_vol=0;
%% Calculation of the actual index

for i=1:nb

    
% Index equation base on average value    
indxs =wv_i(1,i)*((abs(Vnom(i)-Vmean(:,i))./delV(i)).^p_2)';


% Applying considerations to understand the value of the index
i2(i,1)=max(indxs);


rid(i,1) = i2(i,1).^(1/p_2);

if rid(i)<=1
    ridh(i,1)=1;
    i2h(i,1)=1;

else if (rid(i)>=1) && (rid(i)<=5*(dev*100))
     ridh(i,1)=rid(i);
     i2h(i,1)=i2(i);
    else               
     ridh(i,1)= rid(i);
     i2h(i,1)=i2(i);
    end
end

if (rid(i)>1) && (rid(i)<5*(dev*100))
     over_vol=over_vol+1;
     Vol_violation(over_vol,1)=i;
end
end


% Buses with the voltage less/more than the aceptable min/max voltage


if over_vol==0
    Vol_violation=0;
end

i2_tot=sum(i2);

%%%%%%%%%%%% Under/Over Voltage Index %%%%%%%%%%%
V_index = sum(i2h)/(nb);  % The actual index normalized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Ploting results
%  r1 = find(i2>=1);   % index of lines violating limits
%  nvl= size(r1,1);  % number of lines violating limits
%  r2 = (i2).^(1/p_2);
%  nvlsevere=0; 
% 
% if nvl>0
%     figure;plot(t,BVm(:,r1));axis tight
%     title(['Lines violating limits = ', num2str(nvl),'   Index=',num2str(V_index)])
%     xlabel('Time (sec)')
%     ylabel('Voltage')
%     if nvlsevere>0
%         figure;plot(t,BVm(:,r20));axis tight
%         title(['Lines (',num2str(nvlsevere), ') violating limits ',num2str(dev*100),' %'])
%         xlabel('Time (sec)')
%         ylabel('Voltage')
%     end
%     a= find(i2(r1)==max(i2(r1)));
%     mm=r1(1);
% 
%     figure;plot(t,BVm(:,mm));axis tight
%     hold on; plot(t(end-post0:end),BVm(end-post0:end,mm),'linewidth',3)
%     hold on; horline(t,Vmin(mm),'k:');horline(t,Vmax(mm),'k:')
%     hold on; horline(t,mean(Vpost(:,mm)),'r')
%     title(['Line ', num2str(mm),' with maximum violation,  v=',num2str(V_index)])
%     xlabel('Time (sec)')
%     ylabel('Voltage')
% end 
 
%%
clear k;
tt1=simTime-3*deltime;tt2=simTime-2*deltime;tt3=simTime-deltime;tt4=simTime;
V=BVm;
for k=1:nb
count1=0;count2=0;count3=0;
for i=1:length(t)
if(t(i,1)>tt1 && t(i,1)<=tt2)
count1=count1+1;
 V1(count1,k)=V(i,k);
 t1(count1,k)=t(i,1);
  else if (t(i,1)>tt2 && t(i,1)<=tt3)
  count2=count2+1;
  V2(count2,k)=V(i,k);
  t2(count2,k)=t(i,1);
 else if (t(i,1)>tt3 && t(i,1)<=tt4)
 count3=count3+1;
 V3(count3,k)=V(i,k);
 t3(count3,k)=t(i);
     end
      end 
end
end
% Calculating the slopes of the intervals for 'k-th' bus
for i=2:count1
Slope1(i-1,k)=(V1(i,k)-V1(i-1,k))/(t1(i,k)-t1(i-1,k));
if Slope1(i-1,k)== Inf || Slope1(i-1,k)==-Inf 
Slope1(i-1,k)=0;
end
end

for i=2:count2
Slope2(i-1,k)=(V2(i,k)-V2(i-1,k))/(t2(i,k)-t2(i-1,k));
if Slope2(i-1,k)== Inf || Slope2(i-1,k)==-Inf
Slope2(i-1,k)=0;
end
end

for i=2:count3
Slope3(i-1,k)=(V3(i,k)-V3(i-1,k))/(t3(i,k)-t3(i-1,k));
if Slope3(i-1,k)== Inf || Slope3(i-1,k)==-Inf
Slope3(i-1,k)=0;
end
end

% Calculation of mean slope of the interval for 'k-th' bus
 mean_slope(1,k)=mean(Slope1(:,k));
 mean_slope(2,k)=mean(Slope2(:,k));
 mean_slope(3,k)=mean(Slope3(:,k));
 
% Calculation of slope variation between the intervals for 'k-th' bus 
 slope_change1(k,1)=(mean_slope(2,k)-mean_slope(1,k));
 slope_change1(k,2)=(mean_slope(3,k)-mean_slope(2,k));

% Calculation of the difference between the slope variation for 'k-th' bus  
variation(k,1)=(abs(slope_change1(k,2)-slope_change1(k,1)));
end
% Assigning all the buses to an array
 for m=1:nb
 buses(m,1)=m;
 end

 G=[buses variation]; % Matrix with slope variation of bus voltages

% loop to find the lines in which the voltage variation is high
c1=0;   %counting variable
for i=1:nb
if (variation(i)>1)
    c1=c1+1;
    GG(c1,1)=variation(i);
    gg(c1,1)=buses(i);
end
end


% if all voltages reaches steady state then setting g to zero
if c1>0
g=[gg,GG];
else
g=[0,0];
end

% Finding the lines in which reached steady state
SS_bus=0;   % variable to count the lines which reached steady state
for i=1:nb
if ( variation(i,1)<=0.001)
    SS_bus=SS_bus+1;
    steady_state_bus(SS_bus,1)=i;
end

end



end  % function end


 
function   [J]=dynamic_transient(t,t1,delta,M)
%% Integral Square Generator Angle (ISGA) index. 
% Used to evaluate the transient stability of the system, simple approach that looks at the integral
% of the square of the angular deviation from equilibrium (COI).
%
%  [J]=dynamic_transient(t,t1,delta,M)
% 
%  t     - Time vector
%  t1    - [tnom, tend, dt] 
%           tnom - Initial time to be analized (sec) 
%           tend - Final time to be analized (sec)
%             dt - Number of sampling times before the final time, to be analyzed
%  delta - Machine angles vector 
%  M     - Mechanical Starting time, two time interia (H) of the machines (M=2*H)
%
%                      Version 1.1
%


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
 % dt = t1(3);
  
  if  size(M,2)>1
      M=M';
  end

 
%%

nm=size(delta,2);   %machines number
nt   = size(tdel,1);   % simulation length

 delnom  = delta(tnom,:);
 del_uns = [];
 del_uns = find(delta(tend,:)<=-10000); %unstable machines
 
 
 if size(del_uns,2)>0
     nmu=size(del_uns,1);     
     delta0=delta(:,[1:del_uns-1,del_uns+1:end]);
     M0=M([1:del_uns-1,del_uns+1:end],1);
     
     if flag1==1
     figure(401);plot(t,delta(:,del_uns));title(['Stop of generating unit ',num2str(del_uns)])
     figure(400);plot(t,delta0);title(['Machines Angle of the working units (',num2str(nm-nmu),')'])
     end
     
     delta=delta0;
     M=M0;
     nm=nm-size(del_uns,2);
     
 end
    
 delta_a=delta([tnom:tend],:);
 
 
 for i=1:nm
     SM(:,i)=M(i)*delta_a(:,i);
 end
 
num_coa   = sum(SM,2);
den_coa   = sum(M);
delta_coa = num_coa/den_coa;

for i=1:nm
     Jg(:,i)=M(i)*(delta_a(:,i)-delta_coa).^2;
     %Jg(:,i)=(delta(:,i)-delta_coa).^2;
     Jcoa(:,i)=trapz(tdel,Jg(:,i));
end


Mtot = sum(M);
T    = nt;%t(end);
Jtot = sum(Jcoa);
MT   = 1/(Mtot*T);
J1    = MT*Jtot;

J=J1/nm;


%keyboard

function   [J,delta_N,xxx]=dynamic_transient(t,delta,M)
%% Integral Square Generator Angle (ISGA) index. 
% Used to evaluate the transient stability of the system, simple approach that looks at the integral
% of the square of the angular deviation from equilibrium (COI). 
%
%  J =dynamic_transient(t,delta,M)
% 
% INPUTS
%  t     - Time vector
%  delta - Machine angles vector 
%  M     - Mechanical Starting time, two time interia (H) of the machines (M=2*H)
%
% OUTPUTS
%  J - Transien stability index
%
%                      Version 1.2
%


  if  size(M,2)>1
      M=M';
  end

 
%%
nm=size(delta,2);     % machines number
nt   = size(t,1);  % simulation length

post0 = 100;
delta_post = delta(end-post0:end,:);
delta_mean = mean(delta_post);

% Remove generators out of service
xxx=[];k=1;
for i=1:nm
 if delta_mean(1,i)<-1e20
     xxx(k,1)=i;
     delta(:,i)=0;
     delta_mean(1,i)=0;
     M(i)=0;
     k=k+1;
 end
end
clear k
delta_N=delta;


 
 for i=1:nm
     SM(:,i)=M(i)*delta(:,i);
 end
 
num_coa   = sum(SM,2);
den_coa   = sum(M);
delta_coa = num_coa/den_coa;

for i=1:nm
     Jg(:,i)=M(i)*(delta(:,i)-delta_coa).^2;
     Jcoa(:,i)=trapz(t,Jg(:,i));
end


Mtot = sum(M);
T    = nt;
Jtot = sum(Jcoa);
MT   = 1/(Mtot*T);
J1   = MT*Jtot;

J=J1/nm;


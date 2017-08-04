function [Pmax,Vmax,P_pvcurve,V_pvcurve,states]=dynamic_voltagestbility_meanvalues(t,Nb,Nl,BVm,LP1,LQ1,Pmax0,x0,rn,flag1,clr1,clr2,lwid,fnum)
%% Estimate Pmax and Vmax
%
%
%   [Ph,Vh,x]=dynamic_voltagestbility(t,Nb,Nl,V,P,Q,x0,r0,flag1,clr1,clr2,lwid)
%  t  - Time vector
%  Nb - Bus number to be analyzed
%  Nl - Line number to be analyzed
%  V  - Bus Voltages vector
%  P  - Active power flows on lines
%  Q  - Reactive pwoer flows on lines
%  x0 - inital guess to compute [d,E,X]
%  r0 - intial guess to compute [alpha,beta] 
%  flag1 - 0 no plots
%          1 shhow plots 
%               clr1 - color for estimated PV data
%               clr2 - color for actual PV data
%               lwid - line width for actial data

%
%%
    

[lsec,nsec]=size(BVm);

jj=1;nn=lsec;
for i=1:nsec
    V1(jj:nn,1)=BVm(:,i);
    Vk(i,1)=mean(BVm(:,i));
    P1(jj:nn,1)=LP1(:,i);
    Pk(i,1)=mean(LP1(:,i));
    Q1(jj:nn,1)=LQ1(:,i);
    Qk(i,1)=mean(LQ1(:,i));
    t1(jj:nn,1)=t(:,i);
    tk(i,1)=mean(t(:,i));
    ns=size(V1,1);
    jj=ns+1;
    nn=nn+lsec;
end


V = Vk;
P = Pk;
Q = Qk;


%%
global V P Q



x01 = x0(1)*ones(size(P));
x02 = mean(V1);
x0 = [x01; x02; x0(3)];

z0 = lsqnonlin(@EstimateE_and_X,x0);

delta = z0(1:length(P));
E  = z0(length(P)+1);
X  = z0(length(P)+2);


r0 = lsqnonlin(@Estimate_a_and_b,rn);

alpha = r0(1);
beta = r0(2);

states=[delta;E;X;alpha;beta];


%--- Calculating P-V curve ---%
[Pmax,Pmaxindex]=max(abs(P));
Pdirection=sign(P(Pmaxindex));
Pmax=0; Pestim=[]; Vs=[]; Vu=[];
maxPflag=0;


while maxPflag==0
    Qestim=alpha+beta*Pmax;
    Vs2=E^2/2-Qestim*X+sqrt(E^4/4-X^2*Pmax^2-Qestim*X*E^2);
    if imag(Vs2)==0 & real(Vs2)>=0
        Vs=[Vs sqrt(Vs2)];
        Vu=[sqrt(E^2/2-Qestim*X-sqrt(E^4/4-X^2*Pmax^2-Qestim*X*E^2)) Vu];
        Pestim=[Pestim Pmax];
        Pmax=Pmax+.01*Pdirection;
    else
        maxPflag=1;
    end
end

[Pmax,Pmaxindex]=max(abs(Pestim));

%--- Save the Pmax and Vmax values ---%
Pmax = Pestim(Pmaxindex);
Vmax = Vs(Pmaxindex);


%--- Save PV curve data for plots ---%
V_pvcurve = [Vs Vu];
if Pdirection >= 0
    P_pvcurve = [Pestim sort(Pestim,2,'descend')];
else
    P_pvcurve = [Pestim sort(Pestim,2,'ascend')];
end


%% 
mm=exist('clr1');
nn=exist('lwid');
kk=exist('fnum');
if flag1==1

    if mm
           figure;plot(P,V,clr2,'linewidth',lwid)
           hold on;plot(P_pvcurve,V_pvcurve,clr1);
           legend('actual','estimated')
           hold on; verline(Pmax,clr1)
           %hold on; hline(Vmax,clr1)
           hold on; verline(Pmax-(Pmax*Pmax0)/100,clr1)
    else       
            
            
        figure;plot(P,V,'r','linewidth',3)
        hold on;plot(P_pvcurve,V_pvcurve,'b');
        legend('actual','estimated')
        hold on; verline(Pmax,'b')
        %hold on; hline(Vmax,'b')
        hold on; verline(Pmax-(Pmax*Pmax0)/100,'b:')
    end
    
   if kk
       figure(fnum);hold on;plot(P,V,clr2,'linewidth',lwid)
           hold on;plot(P_pvcurve,V_pvcurve,clr1);
           legend('actual','estimated')
           hold on; verline(Pmax,clr1)
           %hold on; hline(Vmax,clr1)
           hold on; verline(Pmax-(Pmax*Pmax0)/100,clr1)
   end
       

end
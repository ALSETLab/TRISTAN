function [Pmax,Vmax,P_pvcurve,V_pvcurve,state,E_post,X_post,PPost,VPost,alp,bet]=dynamic_voltagestability_postcontingency(t,P,V,Q,faultime,deltime,simTime,samples)

mm=size(V);
for m=1:mm(2)
s=5;
P_Post1=zeros(s,m); P_Post2=zeros(s,m); P_Post3=zeros(s,m); P_Post4=zeros(s,m); P_Post5=zeros(s,m); 
V_Post1=zeros(s,m); V_Post2=zeros(s,m); V_Post3=zeros(s,m); V_Post4=zeros(s,m); V_Post5=zeros(s,m); 
Q_Post1=zeros(s,m); Q_Post2=zeros(s,m); Q_Post3=zeros(s,m); Q_Post4=zeros(s,m); Q_Post5=zeros(s,m); 

%Post-contingency counting variables
count1_post=0; count2_post=0; count3_post=0; count4_post=0; count5_post=0;
%%
max=size(t);
loadvar=round((simTime-faultime)/deltime);
for ii=1:max
    if loadvar>=5
         if (t(ii,1)>=faultime && t(ii,1)<(faultime+deltime))
         count1_post=count1_post+1;
         P_Post1(count1_post,m)=P(ii,m);
         V_Post1(count1_post,m)=V(ii,m);
         Q_Post1(count1_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+deltime) && t(ii,1)<(faultime+2*deltime))
         count2_post=count2_post+1;
         P_Post2(count2_post,m)=P(ii,m);
         V_Post2(count2_post,m)=V(ii,m);
         Q_Post2(count2_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+2*deltime) && t(ii,1)<(faultime+3*deltime))
         count3_post=count3_post+1;
         P_Post3(count3_post,m)=P(ii,m);
         V_Post3(count3_post,m)=V(ii,m);
         Q_Post3(count3_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+3*deltime) && t(ii,1)<(faultime+4*deltime))
         count4_post=count4_post+1;
         P_Post4(count4_post,m)=P(ii,m);
         V_Post4(count4_post,m)=V(ii,m);
         Q_Post4(count4_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+4*deltime) && t(ii,1)<(faultime+5*deltime))
         count5_post=count5_post+1;
         P_Post5(count5_post,m)=P(ii,m);
         V_Post5(count5_post,m)=V(ii,m);
         Q_Post5(count5_post,m)=Q(ii,m);
             end
             end
             end
             end
         end
    end
         
    if (loadvar==4)
         if (t(ii,1)>=faultime && t(ii,1)<(faultime+deltime))
         count1_post=count1_post+1;
         P_Post1(count1_post,m)=P(ii,m);
         V_Post1(count1_post,m)=V(ii,m);
         Q_Post1(count1_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+deltime) && t(ii,1)<(faultime+2*deltime))
         count2_post=count2_post+1;
         P_Post2(count2_post,m)=P(ii,m);
         V_Post2(count2_post,m)=V(ii,m);
         Q_Post2(count2_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+2*deltime) && t(ii,1)<(faultime+3*deltime))
         count3_post=count3_post+1;
         P_Post3(count3_post,m)=P(ii,m);
         V_Post3(count3_post,m)=V(ii,m);
         Q_Post3(count3_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+3*deltime) && t(ii,1)<(faultime+4*deltime))
         count4_post=count4_post+1;
         P_Post4(count4_post,m)=P(ii,m);
         V_Post4(count4_post,m)=V(ii,m);
         Q_Post4(count4_post,m)=Q(ii,m);
             end
             end
             end
         end  
    end
    if (loadvar==3)
         if (t(ii,1)>=faultime && t(ii,1)<(faultime+deltime))
         count1_post=count1_post+1;
         P_Post1(count1_post,m)=P(ii,m);
         V_Post1(count1_post,m)=V(ii,m);
         Q_Post1(count1_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+deltime) && t(ii,1)<(faultime+2*deltime))
         count2_post=count2_post+1;
         P_Post2(count2_post,m)=P(ii,m);
         V_Post2(count2_post,m)=V(ii,m);
         Q_Post2(count2_post,m)=Q(ii,m);
         else if (t(ii,1)>=(faultime+2*deltime) && t(ii,1)<(faultime+3*deltime))
         count3_post=count3_post+1;
         P_Post3(count3_post,m)=P(ii,m);
         V_Post3(count3_post,m)=V(ii,m);
         Q_Post3(count3_post,m)=Q(ii,m);
             end
             end
         end         
    end         
end

diff=10^-3;
countPP1=0;countPP2=0;countPP3=0;countPP4=0;countPP5=0;


if loadvar>=5
 kk=count5_post;
else if loadvar==4
 kk=count4_post; 
    else if loadvar==3
            kk=count3_post;
        end
    end
end

for ii=2:kk-1
if loadvar>=5
    if(P_Post1(ii,m)-P_Post1(ii-1,m)<=diff && countPP1<samples)
        countPP1=countPP1+1;
        PPost1(countPP1,m)=P_Post1(ii-1,m);
        QPost1(countPP1,m)=Q_Post1(ii-1,m);
        VPost1(countPP1,m)=V_Post1(ii-1,m);
    end
    
    if(P_Post2(ii,m)-P_Post2(ii-1,m)<=diff && countPP2<samples)
        countPP2=countPP2+1;
        PPost2(countPP2,m)=P_Post2(ii-1,m);
        QPost2(countPP2,m)=Q_Post2(ii-1,m);
        VPost2(countPP2,m)=V_Post2(ii-1,m);
    end
    
    if(P_Post3(ii,m)-P_Post3(ii-1,m)<=diff && countPP3<samples)
        countPP3=countPP3+1;
        PPost3(countPP3,m)=P_Post3(ii-1,m);
        QPost3(countPP3,m)=Q_Post3(ii-1,m);
        VPost3(countPP3,m)=V_Post3(ii-1,m);
    end   
    
    if(P_Post4(ii,m)-P_Post4(ii-1,m)<=diff && countPP4<samples)
        countPP4=countPP4+1;
        PPost4(countPP4,m)=P_Post4(ii-1,m);
        QPost4(countPP4,m)=Q_Post4(ii-1,m);
        VPost4(countPP4,m)=V_Post4(ii-1,m);
    end
    
    if(P_Post5(ii,m)-P_Post5(ii-1,m)<=diff && countPP5<samples)
        countPP5=countPP5+1;
        PPost5(countPP5,m)=P_Post5(ii-1,m);
        QPost5(countPP5,m)=Q_Post5(ii-1,m);
        VPost5(countPP5,m)=V_Post5(ii-1,m);
    end
    
end
if loadvar==4
    if(P_Post1(ii,m)-P_Post1(ii-1,m)<=diff && countPP1<samples)
        countPP1=countPP1+1;
        PPost1(countPP1,m)=P_Post1(ii-1,m);
        QPost1(countPP1,m)=Q_Post1(ii-1,m);
        VPost1(countPP1,m)=V_Post1(ii-1,m);
    end
    
    if(P_Post2(ii,m)-P_Post2(ii-1,m)<=diff && countPP2<samples)
        countPP2=countPP2+1;
        PPost2(countPP2,m)=P_Post2(ii-1,m);
        QPost2(countPP2,m)=Q_Post2(ii-1,m);
        VPost2(countPP2,m)=V_Post2(ii-1,m);
    end
    
    if(P_Post3(ii,m)-P_Post3(ii-1,m)<=diff && countPP3<samples)
        countPP3=countPP3+1;
        PPost3(countPP3,m)=P_Post3(ii-1,m);
        QPost3(countPP3,m)=Q_Post3(ii-1,m);
        VPost3(countPP3,m)=V_Post3(ii-1,m);
    end   
    
    if(P_Post4(ii,m)-P_Post4(ii-1,m)<=diff && countPP4<samples)
        countPP4=countPP4+1;
        PPost4(countPP4,m)=P_Post4(ii-1,m);
        QPost4(countPP4,m)=Q_Post4(ii-1,m);
        VPost4(countPP4,m)=V_Post4(ii-1,m);
    end     
end

if loadvar==3
     if(P_Post1(ii,m)-P_Post1(ii-1,m)<=diff && countPP1<samples)
        countPP1=countPP1+1;
        PPost1(countPP1,m)=P_Post1(ii-1,m);
        QPost1(countPP1,m)=Q_Post1(ii-1,m);
        VPost1(countPP1,m)=V_Post1(ii-1,m);
    end
    
    if(P_Post2(ii,m)-P_Post2(ii-1,m)<=diff && countPP2<samples)
        countPP2=countPP2+1;
        PPost2(countPP2,m)=P_Post2(ii-1,m);
        QPost2(countPP2,m)=Q_Post2(ii-1,m);
        VPost2(countPP2,m)=V_Post2(ii-1,m);
    end
    
    if(P_Post3(ii,m)-P_Post3(ii-1,m)<=diff && countPP3<samples)
        countPP3=countPP3+1;
        PPost3(countPP3,m)=P_Post3(ii-1,m);
        QPost3(countPP3,m)=Q_Post3(ii-1,m);
        VPost3(countPP3,m)=V_Post3(ii-1,m);
    end
end
end
%%
global PP VV QQ alpha beta
if loadvar>=5
P_Post(:,m)=[PPost1(:,m);PPost2(:,m);PPost3(:,m);PPost4(:,m);PPost5(:,m)];
V_Post(:,m)=[VPost1(:,m);VPost2(:,m);VPost3(:,m);VPost4(:,m);VPost5(:,m)];
Q_Post(:,m)=[QPost1(:,m);QPost2(:,m);QPost3(:,m);QPost4(:,m);QPost5(:,m)];
PPost(:,m)=[mean(PPost3(:,m));mean(PPost4(:,m));mean(PPost5(:,m))];
VPost(:,m)=[mean(VPost3(:,m));mean(VPost4(:,m));mean(VPost5(:,m))];
else if loadvar==4
P_Post(:,m)=[PPost1(:,m);PPost2(:,m);PPost3(:,m);PPost4(:,m)];
V_Post(:,m)=[VPost1(:,m);VPost2(:,m);VPost3(:,m);VPost4(:,m)];
Q_Post(:,m)=[QPost1(:,m);QPost2(:,m);QPost3(:,m);QPost4(:,m)];
PPost(:,m)=[mean(PPost2(:,m));mean(PPost3(:,m));mean(PPost4(:,m))];
VPost(:,m)=[mean(VPost2(:,m));mean(VPost3(:,m));mean(VPost4(:,m))];             
    else if loadvar==3
P_Post(:,m)=[PPost1(:,m);PPost2(:,m);PPost3(:,m)];
V_Post(:,m)=[VPost1(:,m);VPost2(:,m);VPost3(:,m)];
Q_Post(:,m)=[QPost1(:,m);QPost2(:,m);QPost3(:,m)];
PPost(:,m)=[mean(PPost1(:,m));mean(PPost2(:,m));mean(PPost3(:,m))];
VPost(:,m)=[mean(VPost1(:,m));mean(VPost2(:,m));mean(VPost3(:,m))];                             
        end
    end
end
%%
PP=P_Post(:,m);
VV=V_Post(:,m);
QQ=Q_Post(:,m);
y03=0.01;                        %guess value for X
%%
y01 = ones(size(P_Post(:,m)));       %Fixing the size for deltas with guess values
y02 = mean(V_Post(:,m));             %guess value for E
y0 = [y01; y02; y03];
z02(:,m) = lsqnonlin(@Estimate,y0);

delta(:,1) = z02(1:length(PP));
E    = z02(length(PP)+1);
X    = z02(length(PP)+2);
E_post(m,1)=E;
X_post(m,1)=X;
%%
r01 = 1;  % Inital guess for alpha
r02 = 1;  % Initial guess for beta
rn=[r01;r02];
r0(:,m) = lsqnonlin(@Estimate_a_b,rn);
alpha = r0(1);
beta = r0(2);
alp(m)=alpha;
bet(m)=beta;
state(:,m)=[delta;E;X;alpha;beta];
%%

Pdirection=1;
Pmax=0; Pestim=[]; Vs=[]; Vu=[];
maxPflag=0;
countP_esti=0;
while maxPflag==0
    Qestim=alpha+beta*Pmax;
    Vs2=E^2/2-Qestim*X+sqrt(E^4/4-X^2*Pmax^2-Qestim*X*E^2);
    if imag(Vs2)==0 && real(Vs2)>=0
        countP_esti=countP_esti+1;
        Vs=[Vs sqrt(Vs2)];
        Vu=[sqrt(E^2/2-Qestim*X-sqrt(E^4/4-X^2*Pmax^2-Qestim*X*E^2)) Vu];
        Pestim=[Pestim Pmax];
        Pmax=Pmax+.01*Pdirection;
    else
        maxPflag=1;
    end
end


count=1;
for ii=2:countP_esti
if(Pestim(ii)>Pestim(count))
    Pmax(m)=Pestim(ii);
    Pmaxindex=ii;
    count=count+1;
else
    count=count+1;
end
end

%--- Save the Pmax and Vmax values ---%
Vmax(m) = Vs(Pmaxindex);


%--- Save PV curve data for plots ---%
V_pvcurve(m,:) = [Vs Vu];
if Pdirection >= 0
    P_pvcurve(m,:) = [Pestim sort(Pestim,2,'descend')];
else
    P_pvcurve(m,:) = [Pestim sort(Pestim,2,'ascend')];
end

  Pmax0=4;
  figure(m);plot(PP,VV,'r','linewidth',3)
  hold on;plot(P_pvcurve,V_pvcurve,'b');
  legend('actual','estimated')
  hold on; verline(Pmax(m),'b')
  hold on; verline(Pmax(m)-(Pmax(m)*Pmax0)/100,'b:')
end  
end

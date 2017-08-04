function [PPremax1,VPremax1,PPremax2,VPremax2,PPre_pvcurve1,VPre_pvcurve1,PPre_pvcurve2,VPre_pvcurve2,states1,states2,E_pre,X_pre,P1,V1]=dynamic_voltagestability_precontingency(t,P,V,Q,P0,Q0,V0,P2,V2,Q2,faultime,alp,bet,samples)
global PP QQ VV
%Pre-contingency loading conditions
mm=size(V);
for m=1:mm(2) 
s=5;
PPre1=zeros(s,m);
VPre1=zeros(s,m);
QPre1=zeros(s,m);
count_pre=0;    %counting variable
max=size(t);

for ii=1:max
         if t(ii,1)<faultime
         count_pre=count_pre+1;
         PPre1(count_pre,m)=P(ii,m);
         VPre1(count_pre,m)=V(ii,m);
         QPre1(count_pre,m)=Q(ii,m);
         end
end
diff=10^-3;
countpre=0;
for ii=2:count_pre
if(PPre1(ii,m)-PPre1(ii-1,m)<=diff && countpre<samples)
countpre=countpre+1;
PPre(countpre,m)=PPre1(ii-1,m);
QPre(countpre,m)=QPre1(ii-1,m);
VPre(countpre,m)=VPre1(ii-1,m);
end
end

P1=mean(PPre(:,m));
Q1=mean(QPre(:,m));
V1=mean(VPre(:,m));

for ii=1:2

P_Pre(:,m)=[P0(1,m);P1;P2(ii,m)];
V_Pre(:,m)=[V0(1,m);V1;V2(ii,m)];
Q_Pre(:,m)=[Q0(1,m);Q1;Q2(ii,m)];
PP=P_Pre(:,m);
VV=V_Pre(:,m);
QQ=Q_Pre(:,m);
x03=0.1;
x01 = ones(size(P_Pre(:,m)));       %Fixing the size for deltas with guess values
x02 = mean(V_Pre(:,m));             %guess value for E
x0 = [x01; x02; x03];
z01(:,m) = lsqnonlin(@Estimate,x0);

delta(:,1) = z01(1:length(PP));
E  = z01(length(PP)+1);
X  = z01(length(PP)+2);
X_pre(ii,m)=X;
E_pre(ii,m)=E;
if ii==1
states1(:,m)=[delta;E;X];
else
states2(:,m)=[delta;E;X];   
end

Pdirection=1;
Pmax=0; Pestim=[]; Vs=[]; Vu=[];
maxPflag=0;
countP_esti=0;
while maxPflag==0
    Qestim=alp(m)+bet(m)*Pmax;
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
for ll=2:countP_esti
if(Pestim(ll)>Pestim(ll-1))
    Pmax=Pestim(ll);
    Pmaxindex=ll;
end
end
%--- Save the Pmax and Vmax values ---%
Vmax = Vs(Pmaxindex);


%--- Save PV curve data for plots ---%
V_pvcurve= [Vs Vu];
if Pdirection >= 0
    P_pvcurve = [Pestim sort(Pestim,2,'descend')];
else
    P_pvcurve = [Pestim sort(Pestim,2,'ascend')];
end

if ii==1
PPremax1(m)=Pmax;
VPremax1(m)=Vmax;
if m==1
PPre_pvcurve1(m,:)=P_pvcurve;
VPre_pvcurve1(m,:)=V_pvcurve;
kk=length(PPre_pvcurve1(m,:));
else if (m>1 && length(P_pvcurve)==kk)     
  PPre_pvcurve1(m,:)=P_pvcurve;
  VPre_pvcurve1(m,:)=V_pvcurve;
  else if (m>1 && length(P_pvcurve)<kk)
      for aa=length(P_Pvcurve):kk
          P_pvcurve(1,aa)=0;
          V_pvcurve(1,aa)=0;
      end
        PPre_pvcurve1(m,:)=P_pvcurve;
        VPre_pvcurve1(m,:)=V_pvcurve;
      else if (m>1 && length(P_pvcurve)>kk)
       for aa=kk:length(P_pvcurve)
       PPre_pvcurve1(m-1,aa)=0;
       VPre_pvcurve1(m-1,aa)=0;
       end
         PPre_pvcurve1(m,:)=P_pvcurve;
         VPre_pvcurve1(m,:)=V_pvcurve;   
          end
      end
  end
end
clr='-*b';
k='b';
w=3;
else
  PPremax2(m)=Pmax;
  VPremax2(m)=Vmax;
  if m==1
  PPre_pvcurve2(m,:)=P_pvcurve;
  VPre_pvcurve2(m,:)=V_pvcurve;
  bb=length(PPre_pvcurve2);
  else if(m>1 && length(P_pvcurve)==bb)
  PPre_pvcurve2(m,:)=P_pvcurve;
  VPre_pvcurve2(m,:)=V_pvcurve;
  else if (m>1 && length(P_pvcurve)<bb)
      for aa=length(P_Pvcurve):bb
          P_pvcurve(1,aa)=0;
          V_pvcurve(1,aa)=0;
      end
        PPre_pvcurve2(m,:)=P_pvcurve;
        VPre_pvcurve2(m,:)=V_pvcurve;
      else if (m>1 && length(P_pvcurve)>bb)
       for aa=bb:length(P_pvcurve)
       PPre_pvcurve2(m-1,aa)=0;
       VPre_pvcurve2(m-1,aa)=0;
       end
         PPre_pvcurve2(m,:)=P_pvcurve;
         VPre_pvcurve2(m,:)=V_pvcurve;   
          end
      end
      end
  end
  clr='-*g';
  k='g';
  w=1.5;
end

  Pmax0=4;hold on;
  plot(PP,VV,clr,'linewidth',w)
   hold on (m);plot(PPre(:,m),VPre(:,m),'linewidth',w);
   hold on (m);plot(P_pvcurve,V_pvcurve,k);
 
 legend('actual','estimated')
  hold on (m); verline(Pmax,'b')
  hold on (m); verline(Pmax-(Pmax*Pmax0)/100,'b:')
end
end
end

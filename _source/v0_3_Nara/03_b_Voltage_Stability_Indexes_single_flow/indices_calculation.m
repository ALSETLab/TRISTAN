function [SBI1,ABI1,GBI1,SBI2,ABI2,GBI2,DataPre1,DataPost1,DataPre2,DataPost2]=indices_calculation(PPremax,VPremax,Pmax0_pre,Pmax0_post,Vmax0_pre,P0,P1,P2,V0,V1,V2,P_pvpost,V_pvpost,PPost,VPost,PPostmax,VPostmax,P_pvpre1,V_pvpre1,P_pvpre2,V_pvpre2)
global Lb
for m=1:length(Lb)
for kk=1:2 
Plim_pre=PPremax(kk,m)-(Pmax0_pre(m)*PPremax(kk,m))/100
Vlim_pre=VPremax(kk,m)+(Vmax0_pre(m)*VPremax(kk,m))/100;

Pmean_pre=[P0(1,m);P1;P2(kk,m)]
Vmean_pre=[V0(1,m);V1;V2(kk,m)];

PLimits_pre=(Plim_pre-Pmean_pre)/Plim_pre
VLimits_pre=(Vmean_pre-Vlim_pre)/Vlim_pre;

% For Pre-Contingency PV curve-2
Pmean_post=PPost(:,m)
Vmean_post=VPost(:,m);

Pmax_post=PPostmax(1,m)
Vmax_post=VPostmax(1,m);

Plim_post=Pmax_post-(Pmax_post*Pmax0_post(1,m))/100
Vlim_post=Vmax_post+(Vmax_post*Vmax0_pre(1,m))/100;
PLimits_post=(Plim_post-Pmean_post)/Plim_post;
VLimits_post=(Vmean_post-Vlim_post)/Vlim_post;
%%
lmax=1.0;
lmin=0.9;

Pnom0 = find(P_pvpost(m,:)>=Pmean_pre(1)*lmin & P_pvpost(m,:)<=Pmean_pre(1)*lmax  &  V_pvpost(m,:)>=Vmax_post);
   if size(Pnom0,2)>=3
            Pnom0 = find(P_pvpost(m,:)>=Pmean_pre(1)*0.99 & P_pvpost(m,:)<=Pmean_pre(1)*1.01 &  V_pvpost(m,:)>=Vmax_post);
   end   
   if size(Pnom0,2)==0
            Pnom0 = find(P_pvpost(m,:)>=Pmean_pre(1)*0.9 & P_pvpost(m,:)<=Pmean_pre(1)*1.1 &  V_pvpost(m,:)>=Vmax_post);
   end
   absVal=abs(Pmean_pre(1)-P_pvpost(m,Pnom0));
   minVal=min(abs(Pmean_pre(1)-P_pvpost(m,Pnom0)));
   selVal=find(absVal==minVal);
   Pnom_A=Pnom0(selVal);
   
   %%
   Pnom1 = find(P_pvpost(m,:)>=Pmean_pre(2)*lmin & P_pvpost(m,:)<=Pmean_pre(2)*lmax  &  V_pvpost(m,:)>=Vmax_post);
   if size(Pnom1,2)>=3
            Pnom1 = find(P_pvpost(m,:)>=Pmean_pre(2)*0.99 & P_pvpost(m,:)<=Pmean_pre(2)*1.01 &  V_pvpost(m,:)>=Vmax_post);
   end   
   if size(Pnom1,2)==0
            Pnom1 = find(P_pvpost(m,:)>=Pmean_pre(2)*0.9 & P_pvpost(m,:)<=Pmean_pre(2)*1.1 &  V_pvpost(m,:)>=Vmax_post);
   end
   absVal1=abs(Pmean_pre(2)-P_pvpost(m,Pnom1));
   minVal1=min(abs(Pmean_pre(2)-P_pvpost(m,Pnom1)));
   selVal1=find(absVal1==minVal1);
   Pnom_B=Pnom1(selVal1); 
   
   %%
    
if Pmax_post>= Pmean_pre(3)  
    Pnom2 = find(P_pvpost(m,:)>=Pmean_pre(3)*lmin & P_pvpost(m,:)<=Pmean_pre(3)*lmax  &  V_pvpost(m,:)>=Vmax_post);
   if size(Pnom2,2)>=3
            Pnom2 = find(P_pvpost(m,:)>=Pmean_pre(3)*0.99 & P_pvpost(m,:)<=Pmean_pre(3)*1.01 &  V_pvpost(m,:)>=Vmax_post);
   end   
   if (size(Pnom2,2)==0 )
           Pnom2 = find(P_pvpost(m,:)>=Pmean_pre(3)*0.9 & P_pvpost(m,:)<=Pmean_pre(3)*1.1 &  V_pvpost(m,:)>=Vmax_post);
   end
   absVal2=abs(Pmean_pre(3)-P_pvpost(m,Pnom2));
   minVal2=min(abs(Pmean_pre(3)-P_pvpost(m,Pnom2)));
   selVal2=find(absVal2==minVal2);
   Pnom_C=Pnom2(selVal2); 
   
   PnomC= P_pvpost(m,Pnom_C)
   VnomC= V_pvpost(m,Pnom_C);
 
else
   
   PnomC=Pmean_pre(3);
   VnomC=Vmax_post;
end
   
 PnomA= P_pvpost(Pnom_A);
 VnomA= V_pvpost(Pnom_A);
 
 PnomB= P_pvpost(Pnom_B);
 VnomB= V_pvpost(Pnom_B);
 
 PnomX=[PnomA;PnomB;PnomC];
 VnomX=[VnomA;VnomB;VnomC];
   
%  Vlim_post1=(VnomA-Vlim_pre)/Vlim_pre;
%  Vlim_post2=(VnomB-Vlim_pre)/Vlim_pre;
%  Vlim_post3=(VnomC-Vlim_pre)/Vlim_pre;
 
 Vlim_post1=(VnomA-Vlim_post)/Vlim_post;
 Vlim_post2=(VnomB-Vlim_post)/Vlim_post;
 Vlim_post3=(VnomC-Vlim_post)/Vlim_post;
 
 Vlim_post0=[Vlim_post1;Vlim_post2;Vlim_post3]
 
 Plim_post1=(Plim_post-PnomA)/Plim_post;
 Plim_post2=(Plim_post-PnomB)/Plim_post;
 Plim_post3=(Plim_post-PnomC)/Plim_post;
 
 Plim_post0=[Plim_post1;Plim_post2;Plim_post3]
 %%

SBI=[  PLimits_pre(1), PLimits_pre(2), PLimits_pre(3),VLimits_pre(1),VLimits_pre(2),VLimits_pre(3);...
        Plim_post0(1) , Plim_post0(2) , Plim_post0(3) ,Vlim_post0(1) ,Vlim_post0(2) ,Vlim_post0(3)]
%%
ABI=[min(SBI(:,1)),min(SBI(:,2)),min(SBI(:,3)),min(SBI(:,4)),min(SBI(:,5)),min(SBI(:,6))]

GBI=[min(ABI(1:3)),min(ABI(4:6))]
%%
if kk==1
if m==1
ll=1;
PP=1;
end
SBI1(ll:ll+1,:)=SBI;
ABI1(PP,:)=ABI;
GBI1(PP,:)=GBI;
DataPre1(:,:,m).limits = [Plim_pre;Vlim_pre];
DataPre1(:,:,m).mean  = [P0 P1 P2(kk);V0 V1 V2(kk)];
DataPre1(:,:,m).deltas = [PLimits_pre';VLimits_pre'];
DataPre1(:,:,m).curve = [P_pvpre1',V_pvpre1'];
DataPost1(:,:,m).limits = [Plim_post;Vlim_post];
DataPost1(:,:,m).mean   = [PnomX';VnomX'];
DataPost1(:,:,m).deltas = [Plim_post0';Vlim_post0'];
DataPost1(:,:,m).curve = [P_pvpost',V_pvpost']; 
ll=ll+2;
PP=PP+1;
else
if m==1    
hh=1;
qq=1;
end
SBI2(hh:hh+1,:)=SBI;
ABI2(qq,:)=ABI;
GBI2(qq,:)=GBI;

DataPre2(:,:,m).limits = [Plim_pre;Vlim_pre];
DataPre2(:,:,m).mean   = [P0 P1 P2(kk);V0 V1 V2(kk)];
DataPre2(:,:,m).deltas = [PLimits_pre';VLimits_pre'];
DataPre2(:,:,m).curve  = [P_pvpre2',V_pvpre2'];
DataPost2(:,:,m).limits = [Plim_post;Vlim_post];
DataPost2(:,:,m).mean   = [PnomX';VnomX'];
DataPost2(:,:,m).deltas = [Plim_post0';Vlim_post0'];
DataPost2(:,:,m).curve = [P_pvpost',V_pvpost']; 
hh=hh+2;
qq=qq+1;
end
end

end

end

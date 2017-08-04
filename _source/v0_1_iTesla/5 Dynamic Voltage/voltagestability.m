function [SBI,ABI,GBI,DataPre,DataPost]=voltagestability(t,V,P,Q,t_pre,t_post,Nb,Nl,Limits,x0,r0)

Pmax0_pre  = Limits(1);  %Maximum Power limit in percent pre-fault
Vmax0_pre  = Limits(2);  %Maximum voltage drop in precent pre-fault 
Pmax0_post = Limits(1);  %Maximum Power limit in percent post-fault

t20=cell2mat(t(1)); t60=cell2mat(t(2)); t80=cell2mat(t(3));
V20=cell2mat(V(1)); V60=cell2mat(V(2)); V80=cell2mat(V(3));
P20=cell2mat(P(1)); P60=cell2mat(P(2)); P80=cell2mat(P(3));
Q20=cell2mat(Q(1)); Q60=cell2mat(Q(2)); Q80=cell2mat(Q(3));

Pre_t(1,:)=cell2mat(t_pre(1));    Pos_t(1,:)=cell2mat(t_post(1));
Pre_t(2,:)=cell2mat(t_pre(2));    Pos_t(2,:)=cell2mat(t_post(2));
Pre_t(3,:)=cell2mat(t_pre(3));    Pos_t(3,:)=cell2mat(t_post(3));


qq=size(Nb,2);kk=1;nn=2;


for jj=1:qq

%%

for j=1:2

    if j==1; pre=1; post=0; end 
    if j==2; pre=0; post=1; end



for i=1:3

    if pre
    %%Pre-fault data   
    t1 = Pre_t(i,:);
    if i==1;  t=t20; end
    if i==2;  t=t60; end    
    if i==3;  t=t80; end 
    end
    
    if post  
    %%Post-fault data
    t1 = Pos_t(i,:);
    if i==1;   t=t20; end
    if i==2;   t=t60; end    
    if i==3;   t=t80; end 
    end
    
lmax=1.0;
lmin=0.9;
  
tnom0 = find(t>=t1(1)*lmin & t<=t1(1)*lmax );
   if size(tnom0,1)>=3
            tnom0 = find(t>=t1(1)*0.999 & t<=t1(1)*1.001 );
   end   
   if size(tnom0,1)==0
            tnom0 = find(t>=t1(1)*0.9 & t<=t1(1)*1.1 );
   end
   
tend0 = find(t>=t1(2)*lmin & t<=t1(2)*lmax);
   if size(tend0,1)>=3
            tend0 = find(t>=t1(2)*0.999 & t<=t1(2)*1.001 );
   end  
   if size(tend0,1)==0
            tend0 = find(t>=t1(2)*0.9 & t<=t1(2)*1.1 );
   end

tnom = tnom0(1);
tend = tend0(end);

tdel0 = t(tnom:tend);


if i==1; tdel=[tdel0]; tkk0=tnom:tend; td=mean(tdel); end
if i==2; tdel=[tdel;tdel0]; tkk1=tnom:tend;  td=[td;mean(tdel)];end
if i==3; tdel=[tdel;tdel0]; tkk2=tnom:tend; td=[td;mean(tdel)]; end
end


a0=size(tkk0,2); a1=size(tkk1,2); a2=size(tkk2,2);

n=min([a0,a1,a2]);


if pre

t0_pre   = [t(tkk0(1:n));t(tkk1(1:n));t(tkk2(1:n))];
Vmag_pre = [ V20(tkk0(1:n),Nb(jj));V60(tkk1(1:n),Nb(jj)); V80(tkk2(1:n),Nb(jj))];
Pij_pre  = [ P20(tkk0(1:n),Nl);P60(tkk1(1:n),Nl); P80(tkk2(1:n),Nl)];
Qij_pre  = [  Q20(tkk0(1:n),Nl);Q60(tkk1(1:n),Nl); Q80(tkk2(1:n),Nl)];


tr_pre    = [ mean(t(tkk0(1:n)));mean(t(tkk1(1:n)));mean(t(tkk2(1:n)))];
Vmean_pre = [ mean(V20(tkk0(1:n),Nb(jj)));mean(V60(tkk1(1:n),Nb(jj))); mean(V80(tkk2(1:n),Nb(jj)))];
Pmean_pre = [ mean(P20(tkk0(1:n),Nl)); mean(P60(tkk1(1:n),Nl)); mean(P80(tkk2(1:n),Nl))];
Qmean_pre = [ mean(Q20(tkk0(1:n),Nl)); mean(Q60(tkk1(1:n),Nl)); mean(Q80(tkk2(1:n),Nl))];

tm_pre=[t(tkk0(1:n)),t(tkk1(1:n)),t(tkk2(1:n))];
Vm_pre=[V20(tkk0(1:n),Nb(jj)),V60(tkk1(1:n),Nb(jj)),V80(tkk2(1:n),Nb(jj))];
Pm_pre=[P20(tkk0(1:n),Nl),P60(tkk1(1:n),Nl),P80(tkk2(1:n),Nl)];
Qm_pre=[Q20(tkk0(1:n),Nl),Q60(tkk1(1:n),Nl),Q80(tkk2(1:n),Nl)];

% x01 = pi/8;             % inital guess for deltas  
% x02 = mean(Vm_pre(:,1));% intial guess for voltage  E
% x03 = 0.03;             % intial guess for reactance X
% 
% r01 = 1;  % Inital guess alpha
% r02 = 1;  % Initial guess beta
% 
% x0 = [x01; x02; x03];
% r0 = [r01,r02];

clear global
[Pmax_pre,Vmax_pre,P_pvpre,V_pvpre,x]=dynamic_voltagestbility_meanvalues(tm_pre,Nb(jj),Nl,Vm_pre,Pm_pre,Qm_pre,Pmax0_pre,x0,r0,0,'b-','r:*',2);

Plim_pre=Pmax_pre-(Pmax0_pre*Pmax_pre)/100;
Vlim_pre=Vmax_pre+(Vmax0_pre*Vmax_pre)/100;
PLimits_pre=(Plim_pre-Pmean_pre)/Plim_pre;
VLimits_pre=(Vmean_pre-Vlim_pre)/Vlim_pre;

DataPre(:,:,jj).limits = [Plim_pre;Vlim_pre];
DataPre(:,:,jj).mean   = [Pmean_pre';Vmean_pre'];
DataPre(:,:,jj).deltas = [PLimits_pre';VLimits_pre'];
DataPre(:,:,jj).curve = [P_pvpre',V_pvpre'];



end

if post

 
t0_post   = [t(tkk0(1:n));t(tkk1(1:n));t(tkk2(1:n))];
Vmag_post = [ V20(tkk0(1:n),Nb(jj));V60(tkk1(1:n),Nb(jj)); V80(tkk2(1:n),Nb(jj))];
Pij_post  = [ P20(tkk0(1:n),Nl);P60(tkk1(1:n),Nl); P80(tkk2(1:n),Nl)];
Qij_post  = [  Q20(tkk0(1:n),Nl);Q60(tkk1(1:n),Nl); Q80(tkk2(1:n),Nl)];


tr_post    = [ mean(t(tkk0(1:n)));mean(t(tkk1(1:n)));mean(t(tkk2(1:n)))];
Vmean_post = [ mean(V20(tkk0(1:n),Nb(jj)));mean(V60(tkk1(1:n),Nb(jj))); mean(V80(tkk2(1:n),Nb(jj)))];
Pmean_post = [ mean(P20(tkk0(1:n),Nl)); mean(P60(tkk1(1:n),Nl)); mean(P80(tkk2(1:n),Nl))];
Qmean_post = [ mean(Q20(tkk0(1:n),Nl)); mean(Q60(tkk1(1:n),Nl)); mean(Q80(tkk2(1:n),Nl))];

tm_post=[t(tkk0(1:n)),t(tkk1(1:n)),t(tkk2(1:n))];
Vm_post=[V20(tkk0(1:n),Nb(jj)),V60(tkk1(1:n),Nb(jj)),V80(tkk2(1:n),Nb(jj))];
Pm_post=[P20(tkk0(1:n),Nl),P60(tkk1(1:n),Nl),P80(tkk2(1:n),Nl)];
Qm_post=[Q20(tkk0(1:n),Nl),Q60(tkk1(1:n),Nl),Q80(tkk2(1:n),Nl)];

% x01 = pi/8;               % inital guess for deltas  
% x02 = mean(Vm_post(:,1)); % intial guess for voltage  E
% x03 = 0.03;               % intial guess for reactance X
% 
% r01 = 1;  %Inital guess alpha
% r02 = 1;  %Initial guess beta
% 
% x0 = [x01; x02; x03];
% r0 = [r01,r02];

clear global
[Pmax_post,Vmax_post,P_pvpost,V_pvpost,x]=dynamic_voltagestbility_meanvalues(tm_post,Nb(jj),Nl,Vm_post,Pm_post,Qm_post,Pmax0_post,x0,r0,0,'b-','r:*',2);

Plim_post=Pmax_post-(Pmax_post*Pmax0_post)/100;
Vlim_post=Vlim_pre;
PLimits_post=(Plim_post-Pmean_post)/Plim_post;
VLimits_post=(Vmean_post-Vlim_post)/Vlim_post;

end


end


%%

lmax=1.0;
lmin=0.9;
  
Pnom0 = find(P_pvpost>=Pmean_pre(1)*lmin & P_pvpost<=Pmean_pre(1)*lmax  &  V_pvpost>=Vmax_post);
   if size(Pnom0,2)>=3
            Pnom0 = find(P_pvpost>=Pmean_pre(1)*0.99 & P_pvpost<=Pmean_pre(1)*1.01 &  V_pvpost>=Vmax_post);
   end   
   if size(Pnom0,2)==0
            Pnom0 = find(P_pvpost>=Pmean_pre(1)*0.9 & P_pvpost<=Pmean_pre(1)*1.1 &  V_pvpost>=Vmax_post);
   end
   absVal=abs(Pmean_pre(1)-P_pvpost(Pnom0));
   minVal=min(abs(Pmean_pre(1)-P_pvpost(Pnom0)));
   selVal=find(absVal==minVal);
   Pnom_A=Pnom0(selVal);
   
Pnom1 = find(P_pvpost>=Pmean_pre(2)*lmin & P_pvpost<=Pmean_pre(2)*lmax  &  V_pvpost>=Vmax_post);
   if size(Pnom1,2)>=3
            Pnom1 = find(P_pvpost>=Pmean_pre(2)*0.99 & P_pvpost<=Pmean_pre(2)*1.01 &  V_pvpost>=Vmax_post);
   end   
   if size(Pnom1,2)==0
            Pnom1 = find(P_pvpost>=Pmean_pre(2)*0.9 & P_pvpost<=Pmean_pre(2)*1.1 &  V_pvpost>=Vmax_post);
   end
   absVal1=abs(Pmean_pre(2)-P_pvpost(Pnom1));
   minVal1=min(abs(Pmean_pre(2)-P_pvpost(Pnom1)));
   selVal1=find(absVal1==minVal1);
   Pnom_B=Pnom1(selVal1);   
   
if Pmax_post>= Pmean_pre(3)
   
    Pnom2 = find(P_pvpost>=Pmean_pre(3)*lmin & P_pvpost<=Pmean_pre(3)*lmax  &  V_pvpost>=Vmax_post);
   if size(Pnom2,2)>=3
            Pnom2 = find(P_pvpost>=Pmean_pre(3)*0.99 & P_pvpost<=Pmean_pre(3)*1.01 &  V_pvpost>=Vmax_post);
   end   
   if (size(Pnom2,2)==0 )
           Pnom2 = find(P_pvpost>=Pmean_pre(3)*0.9 & P_pvpost<=Pmean_pre(3)*1.1 &  V_pvpost>=Vmax_post);
   end
   
   absVal2=abs(Pmean_pre(3)-P_pvpost(Pnom2));
   minVal2=min(abs(Pmean_pre(3)-P_pvpost(Pnom2)));
   selVal2=find(absVal2==minVal2);
   Pnom_C=Pnom2(selVal2);   
   
   PnomC= P_pvpost(Pnom_C);
   VnomC= V_pvpost(Pnom_C);
 
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
 
 Vlim_post1=(VnomA-Vlim_pre)/Vlim_pre;
 Vlim_post2=(VnomB-Vlim_pre)/Vlim_pre;
 Vlim_post3=(VnomC-Vlim_pre)/Vlim_pre;
 
 Vlim_post0=[Vlim_post1;Vlim_post2;Vlim_post3];
 
 Plim_post1=(Plim_post-PnomA)/Plim_post;
 Plim_post2=(Plim_post-PnomB)/Plim_post;
 Plim_post3=(Plim_post-PnomC)/Plim_post;
 
 
  Plim_post0=[Plim_post1;Plim_post2;Plim_post3];

DataPost(:,:,jj).limits = [Plim_post;Vlim_post];
DataPost(:,:,jj).mean   = [PnomX';VnomX'];
DataPost(:,:,jj).deltas = [Plim_post0';Vlim_post0'];
DataPost(:,:,jj).curve = [P_pvpost',V_pvpost'];  
  
%%
SBI([kk,nn],:)=[ PLimits_pre(1), PLimits_pre(2), PLimits_pre(3),VLimits_pre(1),VLimits_pre(2),VLimits_pre(3);...
      Plim_post0(1),Plim_post0(2),Plim_post0(3),Vlim_post0(1),Vlim_post0(2),Vlim_post0(3)];
    kk=kk+2;nn=nn+2;

end
ABI=[min(SBI(:,1)),min(SBI(:,2)),min(SBI(:,3)),min(SBI(:,4)),min(SBI(:,5)),min(SBI(:,6))];

GBI=[min(ABI(1:3)),min(ABI(4:6))];

function [GMI_m,AMI_m,SMI,Poles,Frequencies,Ener]=smallsignal(flagMethod,t,SigAn,t1,Pcent,f,Nf,Damp,flag1);


%% Sorting signals according to energy

nsig = size(SigAn,2);  % number of signals
[ys,tss,ysec,tsec]=energy_sort(t,t1,SigAn,Pcent,flag1);

nsec = size(ysec,2);  % number of signals

if nsec>=100
   y=ysec(:,1:100);
else
y=ysec;
end

Ener.y=ys;
Ener.yh=ysec;
Ener.t=tss;
Ener.th=tsec;

%keyboard

%% Applying ideal filter (FFT) to signals with largest oscillation content

nsig2 = size(y,2); %Number of signals
fmin=f(1);         %Minimum frequency Hz
fmax=f(2);         %Max requency Hz


r=1;
for i=1:nsig2
[f_k,gain_k,]=fft_get(tsec,y(:,i),0.6);

   nfq     = size(f_k,1);
    q      = r+nfq-1;
fqs(r:q,1) = f_k;
fgains(r:q,1)=gain_k;
 fidx(i,1) = nfq; 
   new     = size(fqs,1);
    r      = new+1;
end

Fred_find  = fqs;
nFreq=sum(fidx);
clear fqs

nsig;        % Number of signals analyzed
Fred_find;   % Frequencies bwtween w1 and w2 within the set of signals
fidx;        % Number of frequencies found on each signal  

fss=unique(Fred_find);
Frequencies = fss(find(fss~=0));


%% ERA
if flagMethod==1 

if isempty(Nf)
Nf=size(Frequencies,1);
end



nsig_y=size(y,2);
fr   =  Nf;
n=2*fr+2;          %number of singular values to evaluate
T    =  t(10)-t(9);
[sysd,sv,mc] =  eraiTesla(y,nsig_y,n,T,flag1);

Poles.sys=sysd;
Poles.sv=sv;
Poles.mc=mc;

     nM      = size(sysd.a,1);
     [Aq]    = eig(sysd.a);
     fqmodes = Aq;
     ffreqq  = abs(imag(fqmodes))/2/pi;
     fnfreqq = abs(fqmodes);
     fdampq  = -cos(atan2(imag(fqmodes),real(fqmodes)));
end

%% Prony 

if flagMethod==2 

if isempty(Nf)
Nf=size(Frequencies,1);
end

fr   = Nf;         % Number of modes 
t0   = tsec;       % time vector of order 1xt
yr   = y;          % ringdown matrix of order txNsig, each column is a different signal
Nsig = size(yr,2); % number of signals or channels
n=2*fr+2;          %number of singular values to evaluate

tstart = tsec(1)*ones(1,Nsig);  %starting time for analysis
tend   = tsec(end)*ones(1,Nsig);%ending times for analysis
T      = tsec(end)-tsec(end-1); %sample period for analysis
flag0  = 1;                     %if flag1 = 1, residues are shifted to t = 0.
tstplot= tsec(1);               %starting time for model simulation
tedplot= tsec(end);             %ending time for model simulation
flagN  = flag1;                     %if flag2=1, plot results

[lamda,modelb]=pronyiTesla(t0,yr,n,tstart,tend,T,flag0,tstplot,tsec(end),flagN);

Poles.sys=diag(lamda);

     nM      = size(lamda,1);
     fqmodes = lamda;
     ffreqq  = abs(imag(fqmodes))/2/pi;
     fnfreqq = abs(fqmodes);
     fdampq  = -cos(atan2(imag(fqmodes),real(fqmodes)));
end
     
%% Small Signal Analysis

[frso,idxef] = sort(ffreqq,'ascend');


fqmodes = fqmodes(idxef);
ffreqq  = ffreqq(idxef);
fdampq  = fdampq(idxef);


if flagMethod==1 
fprintf(1,'\n         All Estimated Modes using ERA  ');
end
if flagMethod==2
fprintf(1,'\n         All Estimated Modes using PRONY  ');
end
fprintf(1,'\n    ------------------------------------------ ');
 header = sprintf('\n %15s %17s %10s ',...
        'mode','freq(Hz)','damp (%)');
    disp(header)
    for i=1:nM
        if imag(fqmodes(i))~=0
%            rr=mod(i,2);
%             if rr==1
    tstr = sprintf('%10.4f i%8.4f %10.4f %10.4f',...
     real(fqmodes(i)),imag(fqmodes(i)),ffreqq(i),fdampq(i)*100);
%             else
%                 continue
%             end
        else
            tstr = sprintf('%10.4f       ----- %10.4f %10.4f',...
     real(fqmodes(i)), ffreqq(i),fdampq(i)*100);
        end
     disp(tstr)
    end


ww=find((ffreqq>=fmin) & (ffreqq<=fmax) );

nModes   = size(ww,1);
nMode   = nModes/2;
mods     = fqmodes(ww);
frqmode  = ffreqq(ww);
dampmode = fdampq(ww); 


fprintf(1,'\n              Modes of interest within  ');
fprintf(1,'\n          frequencies: %2.2f and %2.2f Hz ',fmin,fmax);
fprintf(1,'\n    ------------------------------------------ ');
 header = sprintf('\n %15s %17s %10s ',...
        'mode','freq(Hz)','damp (%)');
    disp(header)
    for i=1:nModes
        if imag(mods(i))~=0
%             rr=mod(i,2);
%             if rr==1
     tstr = sprintf('%10.4f i%8.4f %10.4f %10.4f',...
     real(mods(i)),imag(mods(i)),frqmode(i),dampmode(i)*100);
%             else
%                 continue
%             end
        else
            tstr = sprintf('10.4f     ----- %10.4f %10.4f',...
     real(mods(i)), frqmode(i),dampmode(i)*100);
        end
     disp(tstr)
    end
    
dampmode=dampmode*100;

keyboard

th=acos(Damp/100);
ths=pi-th;

thl=acos(dampmode/100);
thls=pi-thl;

jj=1; kk=size(Damp,2);
for i=1:1:nModes%i=1:2:nModes
    for m=1:1:kk        
    SMI(jj,m)=thls(i)-ths(m); % Matrix, Single Mode Indicator(SMI) in rads
    end
    jj=jj+1;
end
 


SMIdeg=SMI*(180/pi);

for m=1:kk
AMI(1,m)=norm(SMI(:,m),2); %Vector, All Mode Indicator (AMI) in rads
AMI_m(1,m)=min(SMI(:,m));
end




  AMIdeg = AMI*(180/pi);
  
  GMI   = norm(AMI,2); % Gain, Global Mode Indicator (GMI) in rads
 GMI_m=min(AMI_m);
  
  



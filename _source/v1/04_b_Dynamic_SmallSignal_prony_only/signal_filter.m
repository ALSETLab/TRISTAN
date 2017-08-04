function [ yh,osc_signal, var_s,ydh ] = signal_filter(sig,t,var_thr,step_min)
%  [Yh,IDX,VAR ]= signal_filter(Y,T,Thr)
%   Filter signal Y based on variance. The signals with threshold lower than Thr are removed. 
%   The time vector T is required, which has variable step size. The output of the function are:
%    Yh- Output signal after filter
%    IDX- indexes of the most relevant signals
%    VAR-Variance of signals Y.
%
% If no oscillation is detected, the function gives empty matrices.



nL  = size(sig,2); % number of signals
ts  = t;           % time vector with variable step
trshold=var_thr;   % select signals with variance higher than var_thr
tsmin= step_min;   % select section of signal with step size grater than step_min 
 
%% Identify disturbance on time-series with variable step

tk=zeros(size(ts));
tq=zeros(size(ts));
tk(2:end)=ts(1:end-1);
tm=ts-tk;
tm(1)=0; %remove first element
tq(3:end)=tm(2:end-1);
h=tm-tq;
n=size(tm,1);
tr=[1:1:n]'; % line improved

%figure;plot(tr,tm)

hmax=max(h);

jj=1; qq=1; ta=zeros(size(tm));

for i=1:n-1,
   ta(jj,:,qq)=tm(i);
   thr0(i,1)=tm(i)-tm(i+1); %
  if tm(i)>hmax
       qq=qq+1;
  end
   jj=jj+1;
end

p=size(ta,3);
ll=size(ta,1);
jj=1;
tw=[];tu=[];
yy=1;
pp=1; ww=0; gg=0;

for i=1:p,
    
for k=1:ll;
       
   if ta(k,:,i)>=0.0001
    tw(jj,1)=ta(k,:,i);
    tw(jj,2)=k;
    jj=jj+1;
   end
   
   if k==ll
    nx=size(tw,1);
    nm(i,1)=size(tw,1)-gg;
    gg=sum(nm);
   end
   
end

 
 for a=yy:nx 
 
     if tw(a,1)<=0.20
         tu(pp,1)=tw(a,1);
         tu(pp,2)=tw(a,2);
     pp=pp+1;
     end
    
   if a==nx
      jo(i,1)=size(tu,1)-ww;
      ww=sum(jo);
   end
 
 end
 
  yy=a+1;
  
  
 
end


%% Signal selection according to variance


time=[]; signal=[]; osc_signal=[];
for i=1:p
 if jo(i)>1e-3
 Var(i,:)=var(sig(tu(sum(jo(1:i-1))+1,2):tu(sum(jo(1:i)),2),:));
 Var_t(i,1)=sum(Var(i,:)); %Variance total
 if Var_t(i,1)>0.01;
     time=ts(tu(sum(jo(1:i-1))+1,2):tu(sum(jo(1:i)),2),:);
     signal=sig(tu(sum(jo(1:i-1))+1,2):tu(sum(jo(1:i)),2),:);   
        
 end   
end
end

if isempty(time)
    osc_signal=[];
    var_s=[];
    yh=[];
    ydh=[];
fprintf('\n Signals not suitable for SSS analysis! \n\n')
else
    j0=1;
    time0=zeros(size(time));
    time0(1)=time(1);
    time0(2:end)=time(1:end-1);
    dt0  =time-time0;
    gg=size(time,1);
    timex=[]; signalx=[];
    for k=1:gg
        if dt0(k)>=tsmin
            dt0x(j0,1)=dt0(k);
            timex(j0,1)=time(k);
            signalx(j0,:)=signal(k,:);
            j0=j0+1;
    end
    end

    if isempty(signalx)
    else
        for k=1:nL,
            signalx0(:,k)=signalx(:,k)-mean(signalx(:,k));
        end
        ydh=[timex,signalx0];
     sig_var=var(signalx0)'; % variance
     sig_std=std(signalx0)'; % standard deviation    
     [ var_sort,s_idx]=sort(sig_var,'descend');
     [ std_sort,std_idx]=sort(sig_std,'descend');
    
     nL0=size(var_sort,1);
     vA=max(var_sort);
     var_sort0=var_sort/vA;
     
     for k=1:nL
        if var_sort0(k)>= trshold;
            osc_signal(k,1)=s_idx(k);
            var_s(k,1)=var_sort(k);
        end
     end
    end
    
    
if isempty(osc_signal)
    osc_signal=[];
    var_s=[];
    yh=[];
    ydh=[];
    fprintf('\n Signals no suitable for SSS analysis! \n\n')
else
  
    
    
 yh=sig(:,osc_signal);
 
 
%  figure;plot(ts,sig(:,osc_signal));
%  title(['Sigansl witht oscillation higher than the threshold (',num2str(trshold),')'])
%  xlabel('Time (sec)')
end
end


function [sysd,sv,mc] =  eraiTesla(y,r,nm,ts,flag1)
% Simple ERA system identification routine 
%  [sysd,sv,mc] = simpera(y,r,nm,ts)
% Outputs:
%  sysd - identified SS discrete LTI system object
%  sv   - singular values (of Hankel Matrix)
%  mc   - modal content = [mode Nfreq(rad/s) freq(Hz) damp residue energy]
% Inputs:
%  y    - impulse response [y1 y2 ... yr], 
%         yi = m impulse responses from ith input
%  r    - number of inputs for MIMO system
%  nm   - number of singular values to evaluate 
%         (note: order of system will be nm*m)
%  ts   - sample time
%


 NChan = size(y,2);
[N,rm] = size(y);     % N = length of data
m = fix(rm/r);        % m = number of outputs
if (m*r ~= rm)
   error('error in number of inputs')
end
t = 0:ts:(N-1)*ts;    % time vector

chansel =1;

% if (m*r ~= 1)   % if it isn't a SISO system
%    for i=1:r   
% 	   for j=1:m
%          disp(['Channel ' num2str( (i-1)*m+j ) '): input ' ...
%                num2str(i) ', output ' num2str(j) '.' ])
%    	end
%    end
%    chansel = Nchannels;
% end

chansel = [1:NChan];

% if flag1==1
 %figstart = figure;
% plot(t,y(:,chansel))
% title('inputed impulse response')
% end


markov = zeros(m,N*r);
for i=1:r         % turn the pulse into Markov Parameters
   mindex = i:r:N*r;
   markov(:,mindex) = y(:,(i-1)*m+1:i*m)';
end

if mod(m/2,2)
   q=(m/2)+0.5;
else
    q=m/2;
end


%or_size=nm*m; %original size
or_size=nm*q; %my size


disp(['Hankel size is ' num2str(or_size) ' by ' num2str((N-nm-1)*r) '.'])

d = markov(1:m,1:r);   % form the D matrix
H0 = zeros(or_size,(N-nm-1)*r);
H1 = zeros(or_size,(N-nm-1)*r);
for i=1:nm             % form the Hankel matrices
   H0((i-1)*m+1:i*m,:) = markov(:,i*r+1:(N-nm-1+i)*r);
   H1((i-1)*m+1:i*m,:) = markov(:,(i+1)*r+1:(N-nm+i)*r);
end

[R,sg,S] = svd(H0);



sgn = sg(:,1:or_size);  %original line
sv = diag(sgn);

% if flag1==1
% figure
% semilogy(sv,'r*'),grid
% title('singular values')
% end



 disp(['The max singular value is ' num2str(sv(1)) ...
       ', and the min is ' num2str(sv(or_size))])
nmt = 0;
for i=1:or_size
   if ((sv(i)/1e-6)>sv(1))
      nmt = nmt+1;
   end
end
nm1 = ceil(nmt/m);   % round up to an even value, new size = nm1*m

% disp(['I suggest using ' num2str(nm1*m) ...
%       ' singular values, then max is ' num2str(sv(1)) ...
%       ' and min is ' num2str(sv(nm1*m))])

%   svflag = input('Choose (1) to agree or (2) to use all: ');

   svflag = 1;



if svflag == 1
	H0n = zeros(nm1*m,(N-nm1-1)*r);
	H1n = zeros(nm1*m,(N-nm1-1)*r);
	for i=1:nm1             % re-form the smaller Hankel matrices
   	H0n((i-1)*m+1:i*m,:) = markov(:,i*r+1:(N-nm1-1+i)*r);
   	H1n((i-1)*m+1:i*m,:) = markov(:,(i+1)*r+1:(N-nm1+i)*r);
	end
	[R,sg,S] = svd(H0n);
else
   nm1 = nm;
   H1n = H1;
end

Rn = R(:,1:nm1*m);   %these are the first nm1 colums of the SVD
sgn = sg(1:nm1*m,1:nm1*m);
Sn = S(:,1:nm1*m);
 
Er = [eye(m) zeros(m,m*(nm1-1))]';     % find the large system matrices
Es = [eye(r) zeros(r,r*(N-nm1-2))]';

Af = sgn^(-0.5)*Rn'*H1n*Sn*sgn^(-0.5);
Bf = sgn^(0.5)*Sn'*Es;
Cf = Er'*Rn*sgn^(0.5);

sysf = ss(Af,Bf,Cf,d,ts);



yf = impulse(sysf,t);
yfp = fiximp(yf,r,m);



% if flag1==1
% figure(figstart)
% plot(t,y(:,chansel),t,yfp(:,chansel),':')
% title('overlaid with large system impulse response')
% end


sysfc = d2c(sysf,'zoh');
[vf,Am] = eig(sysfc.a);
Bm = inv(vf)*sysfc.b;
Cm = sysfc.c*vf;
fmodes = diag(Am);
nm2 = length(fmodes);  % this is in case d2c adds some poles
ffreq = abs(imag(fmodes))/2/pi;
fnfreq = abs(fmodes);
fdamp = -cos(atan2(imag(fmodes),real(fmodes)));

fresidue = zeros(nm2,1);
for i=1:r    % sum over all the inputs and outputs
   for j=1:m
      fresidue = fresidue + abs(Bm(:,i).*Cm(j,:)');
   end
end
fresidue = fresidue/m/r;  % to average it

fenergy = zeros(nm2,1);
for i=1:nm2
   if imag(fmodes(i))==0  % if its a real mode
      fenergy(i) = ts*sum((exp(fmodes(i)*t)).^2);
   else
      Bt = sqrt(1-fdamp(i)^2);
      wt = abs(fmodes(i));  % natural frequency in rad/s
      zt = fdamp(i);   % zeta value
      fenergy(i) = ts*sum((exp(-wt*zt*t)/Bt.*sin(Bt*wt*t+atan2(Bt,zt))).^2);
   end
end




% pre-sorting by frequency
[xx jsort] =sort(ffreq);
while (xx(1)==0)  % for all the real poles
   jsort = [jsort(2:nm2); jsort(1)];  % shift them to the back
   xx = [xx(2:nm2); xx(1)];
end
fmodes = fmodes(jsort);
ffreq = ffreq(jsort);
fnfreq = fnfreq(jsort);
fdamp = fdamp(jsort); 
fresidue = fresidue(jsort);
fenergy = fenergy(jsort);
Bm = Bm(jsort,:); Cm = Cm(:,jsort);

% header = sprintf('\n %6s %14s %10s %6s %13s %10s','mode',...
%    'Nfreq(rad/s)','freq(Hz)','damp','residue','energy');
% disp(header)
% for i=1:nm2
%    tstr = sprintf('%6d %14.4f %10.4f %8.4f %11.4f %11.4f',...
%       i,fnfreq(i),ffreq(i),fdamp(i),fresidue(i),fenergy(i));
%    disp(tstr)
% end

% sortchoice = input('sort by (1)residue or (2)energy: ');
sortchoice = 1;




% sorting
if sortchoice == 1
   [xx jsort] =sort(fresidue);
else
   [xx jsort] =sort(fenergy);
end
jsort = flipud(jsort);
fmodes = fmodes(jsort);
ffreq = ffreq(jsort);
fnfreq = fnfreq(jsort);
fdamp = fdamp(jsort); 
fresidue = fresidue(jsort);
fenergy = fenergy(jsort);
mc = [fmodes fnfreq ffreq fdamp fresidue fenergy];
Bm = Bm(jsort,:); Cm = Cm(:,jsort);

% header = sprintf('\n %6s %14s %10s %6s %13s %10s','mode',...
%    'Nfreq(rad/s)','freq(Hz)','damp','residue','energy');
% disp(header)
% for i=1:nm2
%    tstr = sprintf('%6d %14.4f %10.4f %8.4f %11.4f %11.4f',...
%       i,fnfreq(i),ffreq(i),fdamp(i),fresidue(i),fenergy(i));
%    disp(tstr)
% end

smodes = 1;
selcount = 1;
jay = sqrt(-1);  % ensure imaginary variable
	

while (smodes ~=0)
   
%    smodes = input('select the modes you want as a vector: ');
       
        %npp =2*nm+2; %original line
        npp =nm;  % my line
     smodes = 1:npp;
   
        
     %Check to only have pairs of complex modes
      rr=find(abs(imag(fmodes(smodes)))~=0);
        nrr=size(rr,1);
        if mod(nrr,2) == 0
      smodes = 1:npp;
        else
            if nrr>2
            smodes = 1:npp-1;
            npp=npp-1;
            else
                smodes = 1:npp+1;
                npp=npp+1;
            end
       end 
        
 
    
   
%    if smodes==0
%       break   % break out of this loop
%    end
   
	As = diag(fmodes(smodes));
	Bs = Bm(smodes,:);
	Cs = Cm(:,smodes);

	T=[];   % setup an empty transform matrix
	ns = length(As);
	ii = 1;
	while ii<ns+1
	   if imag(As(ii,ii))==0  % if this is a real mode
         T = [T,zeros(length(T),1);zeros(1,length(T)) 1];  
                     % add a one on diag
	      ii = ii+1;  % advance for while loop
	   elseif ii==ns    % if you try to select a single complex mode
         error('you need to select complex modes in pairs')
      elseif real(As(ii,ii))~=real(As(ii+1,ii+1))  
                     % if you select an unmatched complex mode
         error('you need to select complex modes in pairs')
  		else
     		T = [T,zeros(length(T),2); zeros(2,length(T)) [1 1;jay -jay] ];
	      ii = ii+2;  % skip the conjugate mode
	   end
	end

	Ad = T*As*inv(T);  % do the transformation to get into modal form
	Bd = real(T*Bs);
	Cd = real(Cs*inv(T));
	sysd = ss(Ad,Bd,Cd,d);
	sysdd = c2d(sysd,ts,'zoh');

	yd = impulse(sysdd,t);
	ydp = fiximp(yd,r,m);
	
    if flag1==1
%     figure(figstart+1+selcount)
	%plot(t,y(:,chansel),t,yfp(:,chansel),':',t,ydp(:,chansel),'*')%original
    figure;plot(t,y(:,chansel),t,ydp(:,chansel),'*')
   title(['selected system impulse response, selection number ' ...
        num2str(selcount) ])
    end
    
    
     [Aq] = eig(sysd.a);
     fqmodes = Aq;
     ffreqq = abs(imag(fqmodes))/2/pi;
     fnfreqq = abs(fqmodes);
     fdampq = -cos(atan2(imag(fqmodes),real(fqmodes)));
   
    
    
%     header = sprintf('\n %11s %18s %6s ',...
%         'mode','freq(Hz)','damp');
%     disp(header)
%     for i=1:npp
%         if imag(fqmodes(i))~=0
%     tstr = sprintf('%8.4f i%8.4f %10.4f %8.4f',...
%      real(fqmodes(i)),imag(fqmodes(i)),ffreqq(i),fdampq(i));
%         else
%             tstr = sprintf('%8.4f     ----- %10.4f %8.4f',...
%      real(fqmodes(i)), ffreqq(i),fdampq(i));
%         end
%      disp(tstr)
%     end
    
    
    
    

     
%     disp([' The results of your selection is shown in figure ' ...
%          num2str(figstart+1+selcount) ...
%          ', you can select more modes or enter 0 to quit.' ])
%    selcount = selcount + 1;
 
break   % break out of this loop

end     % end selection while loop

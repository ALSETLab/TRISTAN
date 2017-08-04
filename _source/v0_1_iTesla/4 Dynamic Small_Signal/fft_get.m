function [fsend,magsend]=fft_get(t,x,alpha,flag1);

t=t';
x=x';

n=size(t,2);         % Number of points 

    if mod(n,2) == 0
  n=n;
else
  n=n-1;
    end 



%% Applying FFT to signal "x" and plotting real and imaginary part
spectrum=fft(x);


%keyboard

%%
%Verifying that first and middle components are real
 DCfisrst=spectrum(1:2);
 indNyq=n/2+1;
 NyquistMiddle=spectrum(indNyq-2:indNyq+2);
 
 
 %Fold the spectrum as algorithtm in Figure 12 of document 
 spectrum = spectrum(1:indNyq);       % truncate
 spectrum = spectrum/n;               % scale
 spectrum(2:end) = 2*spectrum(2:end); % compensate for truncating negative freqs
 magnitude = abs(spectrum);           % abs of complex numbers
 
 dt = t(3)-t(2);
 df = (1/dt)/n; % frequency resolution
 
 freqAx=[0:indNyq-1]'*df;    % values of frequency axis
 
 [gains,indxs]=sort(magnitude','descend');
 gains_norm=gains/gains(1);
  kk=find(gains_norm>=alpha);
   mm=size(kk,1);
 
  gsend=indxs(1:mm);
  fsend= freqAx(gsend);
  magsend=magnitude(gsend)';
  
  %keyboard
%% Various Plots
mm=exist('flag1');

if mm


% Plotting signal "x" and its first 20 points
% figure
% subplot(2,1,1)
% plot(t,x);axis tight;grid on
% ylabel('Units')
% xlabel('Time[s]')
% title('Signal')
% subplot(2,1,2)
% plot(t(1:20),x(1:20),'-o');axis tight;grid on
% ylabel('Units')
% xlabel('Time[s]')
 
% figure
% subplot(2,1,1)
% plot(real(spectrum));axis tight;grid on
% ylabel('Real')
% title('Raw FFT spectrum')
% subplot(2,1,2)
% plot(imag(spectrum));axis tight;grid on
% ylabel('Imaginary')
% xlabel('Index number') 

 figure
 plot(freqAx,magnitude);axis tight;grid on
 xlabel('Frequency[Hz]')
 ylabel('Magnitude[units]')
 [magMax,indMax]=max(magnitude);
 freqMax=freqAx(indMax);
 title(sprintf('Max magnitude=%f @ %fHz, \\Deltaf=%fHz',magMax,freqMax,df))
 


%  figure
%  handlestem=stem(freqAx,magnitude);
%  set(handlestem(1),'Marker','none');grid on
%  xlabel('Frequency [Hz]')
%  ylabel('Magnitude [Hz]')
%  title(sprintf('Max magnitude=%f @ %fHz, \\Deltaf=%fHz',magMax,freqMax,df))
%  
% figure
%  fmax=find(freqAx==freqMax);
%  handlestem=stem(freqAx,magnitude);
%  axis([freqAx(fmax-5) freqAx(fmax+5) 0 max(magnitude)])
%  set(handlestem(1),'Marker','none');grid on
%  xlabel('Frequency [Hz]')
%  ylabel('Magnitude [Hz]')
%  title(sprintf('Max magnitude=%f @ %fHz, \\Deltaf=%fHz',magMax,freqMax,df))

end
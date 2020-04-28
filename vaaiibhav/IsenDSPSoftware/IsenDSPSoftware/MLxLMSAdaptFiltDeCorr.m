function MLxLMSAdaptFiltDeCorr(k,Mu,freq,NoTaps,DeCorrDelay)
global   LMSDeCorDataVec Err FiltOut
% function ML_LMSAdaptFiltDeCorr(k,Mu,freq,NoTaps,DeCorrDelay)
%
% k is an amount of white noise to add to a test
% signal which consists of a sinusoid of frequency freq;
%
% Mu is the LMS update term weight;
%
% NoTaps is the number of taps to use in the LMS adaptive filter;
%
% DeCorrDelay is the number of samples to delay the filter input relative to the Plant or channel
% delay.
%
% A typical call:
%
% MLxLMSAdaptFiltDeCorr(0.1,1,100,10,2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
ftsz1 = 14;
ftsz2 = 13;
SR = 1024;
lenFFT = 2*SR;

LMSDeCorDataVec = [];
Err = [];
FiltOut = [];

k = real(k);
Mu = real(Mu);
freq = real(freq);
NoTaps = floor(abs(NoTaps));

if NoTaps > 100
   NoTaps = 100;
end

DeCorrDelay = floor(abs(DeCorrDelay));
if DeCorrDelay > 100
   DeCorrDelay = 100;
end

t = 0:1/SR:1-1/SR;
LMSDeCorDataVec = k*randn(1,length(t))+ 0.5*sin(2*pi*t*1*freq)+ 0.8*cos(2*pi*t*3*freq);
input =  LMSDeCorDataVec;
 TapWt(1:NoTaps,1:length(LMSDeCorDataVec)) = zeros;
 FiltOut(1,1:length(LMSDeCorDataVec)) = zeros;
 Err(1,1:length(LMSDeCorDataVec)) = zeros;
 n = 1:1:NoTaps;
 
for CurDtaPtr = NoTaps+1:1:length(LMSDeCorDataVec)-DeCorrDelay

 FiltOut(1,CurDtaPtr) = LMSDeCorDataVec(CurDtaPtr-n)*TapWt(n,CurDtaPtr);
 Err(1,CurDtaPtr)= LMSDeCorDataVec(1,CurDtaPtr + DeCorrDelay-1)-FiltOut(1,CurDtaPtr); 
LIMAdjust = (std(LMSDeCorDataVec(1,CurDtaPtr-(NoTaps):CurDtaPtr)))^2; 
 TapWt(n,CurDtaPtr+1)= TapWt(n,CurDtaPtr)+ ...
    2*Mu*Err(1,CurDtaPtr)/(50*LIMAdjust)...
    *LMSDeCorDataVec(1,CurDtaPtr-n)';
end


scnsize = get(0,'ScreenSize');
pos1 = [0.005*scnsize(3), 0.06*scnsize(4), 0.99*scnsize(3), 0.86*scnsize(4)];
pos2 = [0.5*scnsize(3), 0.06*scnsize(4), 0.44*scnsize(3), 0.86*scnsize(4)];
figure(17)
set(17,'Name','LMS Adaptive Filter w/Decorrelating Delay',...
   'Color',[1 1 1],'Position',pos1,'Numbertitle','off','paperpositionmode','auto');
rtpt = scnsize(3)-65;

subplot(2,2,1)
plot(LMSDeCorDataVec(1,1:SR/2),'k')
set(gca,'fontsize',ftsz2)
zzc = max(abs(LMSDeCorDataVec));
ylabel(['Amplitude'],'fontsize',ftsz1)
xlabel(['(a) Sample, Input Data'],'fontsize',ftsz1)
axis([0 SR/2 -zzc zzc])

subplot(2,2,3)
plot(Err(1,1:length(t)),'k')
set(gca,'fontsize',ftsz2)
vcx = max(abs(Err));
ylabel(['Amplitude'],'fontsize',ftsz1)
xlabel(['(c) Sample, Output/Error'],'fontsize',ftsz1)
axis([0 inf -zzc zzc])

ScaleFac = max(abs(LMSDeCorDataVec));
LMSDeCorDataVec = LMSDeCorDataVec/ScaleFac;

uicontrol('Style','push','Callback','InLMSDeCorPlay',...
   'Position',[10 260 40 20],'String','Play','backgroundcolor',[1 1 1]);

uicontrol('Style','push','Callback','OutLMSDeCorPlay',...
   'Position',[10 10 40 20],'String','Play','backgroundcolor',[1 1 1]);

%WinLMSDeCorDataVec(1,1:256) = (boxcar(256)').*LMSDeCorDataVec(1,1:256);
lenFFT = 1024;

subplot(2,2,2)

xfft = 20*log10(abs(fft(input+(10^-6),lenFFT)));
zx = max(xfft);
plot([0:1:lenFFT/2-1]/(lenFFT/2),xfft(1,1:lenFFT/2),'k')
set(gca,'fontsize',ftsz2)
axis([0, 1, -inf, zx+5])
xlabel(['(b) Normalized Frequency'],'fontsize',ftsz1 )
ylabel(['Magnitude, dB'],'fontsize',ftsz1)

subplot(2,2,4)

WinErr = Err; 
theans = 20*log10(abs(fft(WinErr+(10^-6),lenFFT)));
plotlim = max(abs(theans));
set(gca,'fontsize',ftsz2)
plot([0:1:lenFFT/2-1]/(lenFFT/2),theans(1,1:lenFFT/2),'k')
axis([0, 1, -inf, max([1.2*plotlim,  zx+5])])
ylabel(['Magnitude, dB'],'fontsize',ftsz1)
xlabel(['(d) Normalized Frequency'],'fontsize',ftsz1 )

figure(88)
clf
set(88,'color',[1,1,1],'name','Enhanced Periodic Signal Output')
uicontrol('Style','push','Callback','FiltOutPlay',...
   'Position',[10 260 40 20],'String','Play','backgroundcolor',[1 1 1]);

lenFFT = 2048;
subplot(211)
plot(FiltOut,'k')
set(gca,'fontsize',ftsz2)
ylabel(['Amplitude'],'fontsize',ftsz1)
xlabel(['(a) Sample, LMS Filter Output'],'fontsize',ftsz1)
axis([0, length(FiltOut)/2, -inf,inf])

subplot(212)
yfr = 20*log10(abs(fft(FiltOut+(10^-6),lenFFT)));
plotlim = max(abs(yfr));
plot([0:1:lenFFT/2-1]/(lenFFT/2),yfr(1,1:lenFFT/2),'k')
set(gca,'fontsize',ftsz2)
ylabel(['Magnitude, dB'],'fontsize',ftsz1)
xlabel(['(b) Normalized Frequency'],'fontsize',ftsz1)
axis([0, 1, -inf, max(yfr)+5])

FiltOut = (0.99/max(abs(FiltOut)))*FiltOut;






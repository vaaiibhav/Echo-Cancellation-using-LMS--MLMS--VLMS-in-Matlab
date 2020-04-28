function LVxLMSAdaptFiltDeCorr(k,Mu,freq,NoTaps,DeCorrDelay)
% function LVxLMSAdaptFiltDeCorr(k,Mu,freq,NoTaps,DeCorrDelay)
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
% LVxLMSAdaptFiltDeCorr(0.1,0.75,200,20,5)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

global LMSDeCorDataVec
global Err FiltOut

LMSDeCorDataVec = [];
Err = [];

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

LMSDeCorDataVec = [];
Err = [];
FiltOut = [];

SR = 2*1024;
t = 0:1/SR:1-1/SR;

LMSDeCorDataVec = 0.5*sin(2*pi*t*0.97*freq)+ 0.5*cos(2*pi*t*3.02*freq)+k*randn(1,length(t));
input =  LMSDeCorDataVec;
 TapWt(1:NoTaps,1:length(LMSDeCorDataVec)) = 0;
 FiltOut(1,1:length(LMSDeCorDataVec)) = 0;
 Err(1,1:length(LMSDeCorDataVec)) = 0;
 n = 1:1:NoTaps;
 
for CurDtaPtr = NoTaps+1:1:length(LMSDeCorDataVec)-DeCorrDelay

 FiltOut(1,CurDtaPtr) = LMSDeCorDataVec(CurDtaPtr-n)*TapWt(n,CurDtaPtr);
 Err(1,CurDtaPtr)= LMSDeCorDataVec(1,CurDtaPtr + DeCorrDelay-1)-FiltOut(1,CurDtaPtr); 
LIMAdjust = (std(LMSDeCorDataVec(1,CurDtaPtr-(NoTaps):CurDtaPtr)))^2; 
 TapWt(n,CurDtaPtr+1)= TapWt(n,CurDtaPtr)+ ...
    2*Mu*Err(1,CurDtaPtr)/(50*LIMAdjust)...
    *LMSDeCorDataVec(1,CurDtaPtr-n)';
end

figure(17)

subplot(2,2,1)
plot(LMSDeCorDataVec(1,1:SR/2),'k')
zzc = max(abs(LMSDeCorDataVec));
ylabel(['Amplitude'])
xlabel(['(a) Sample, Input Data'])
axis([0, SR/2, -zzc, zzc])

subplot(2,2,3)
plot(Err(1,1:length(t)),'k')
vcx = max(abs(Err));
ylabel(['Amplitude'])
xlabel(['(c) Sample, Output/Error'])
axis([0, inf, -zzc, zzc])

ScaleFac = max(abs(LMSDeCorDataVec));
LMSDeCorDataVec = LMSDeCorDataVec/ScaleFac;

WinLMSDeCorDataVec(1,1:SR/4) = LMSDeCorDataVec(1,1:SR/4);
lenFFT = 1024;

subplot(222)

xfft = 20*log10(abs(fft(input+(10^-6),lenFFT)));
zx = max(xfft);
plot([0:1:lenFFT/2-1]/(lenFFT/2),xfft(1,1:lenFFT/2),'k')
axis([0, 1, -inf, zx+5])
xlabel(['(b) Normalized Frequency'])
ylabel(['Magnitude, dB'])

subplot(224)

WinErr = Err; 
theans = 20*log10(abs(fft(WinErr+(10^-6),lenFFT)));
plotlim = max(abs(theans));
plot([0:1:lenFFT/2-1]/(lenFFT/2),theans(1,1:lenFFT/2),'k')
axis([0, 1, -inf, max([1.2*plotlim,  zx+5])])
ylabel(['Magnitude, dB'])
xlabel(['(d) Normalized Frequency'] )

Comment = 'global sound variable names are LMSDeCorDataVec, Err, and FiltOut'
Comment = 'global sound variable sample rate is 8000'
sound( [LMSDeCorDataVec,LMSDeCorDataVec],8000)

figure(88)
clf

lenFFT = 2048;
subplot(211)
plot(FiltOut,'k')
ylabel(['Amplitude'])
xlabel(['(a) Sample, LMS Filter Output'])
axis([0, length(FiltOut)/2, -inf,inf])

subplot(212)
yfr = 20*log10(abs(fft(FiltOut+(10^-6),lenFFT)));
plotlim = max(abs(yfr));
plot([0:1:lenFFT/2-1]/(lenFFT/2),yfr(1,1:lenFFT/2),'k')
ylabel(['Magnitude, dB'])
xlabel(['(b) Normalized Frequency'])
axis([0, 1, -inf, max(yfr)+5])

pause(1.5)

sound([Err,Err],8000)

pause(1.5)

sound([FiltOut,FiltOut],8000)


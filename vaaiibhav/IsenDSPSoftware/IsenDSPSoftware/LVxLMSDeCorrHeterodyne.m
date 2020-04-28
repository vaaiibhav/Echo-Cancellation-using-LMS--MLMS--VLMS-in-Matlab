function LVxLMSDeCorrHeterodyne(strWavFile,A,Mu,Freq,NoTaps,DeCorrDelay)
% function LVxLMSDeCorrHeterodyne(strWavFile,A,Mu,Freq,NoTaps,DeCorrDelay)
% strWavFile is the audio file to open and use as the test signal
% A is the amplitude of an interfering tone (sinusoid) to be added to
% the audio file whose file name or path name (including file name and
% extension) is specified as strWavFile.
% Mu is the LMS update term weight;
% NoTaps is the number of taps to use in the LMS adaptive filter;
% DeCorrDelay is the number of samples to delay the filter input relative to the Plant or channel
% delay.
% The following plots are created: Frequency response (DTFT) of the adaptive filter when converged 
% (magnitude v. frequency in Hz), the test signal (amplitude v. sample),  the filtered test signal (amplitude 
% v. sample), the DTFT of the test signal (magnitude v. frequency in Hz, and the DTFT (magnitude v. 
% frequency in Hz) of the last 50 % of the filtered test signal.
% Test calls:
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.5,0.08,250,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.5,0.08,500,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.5,0.08,1200,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.5,0.01,75,45,10)
%
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',1.5,0.08,250,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',1.5,0.08,500,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',1.5,0.08,1200,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',1.5,0.01,75,45,10)
%
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',4.5,0.08,250,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',4.5,0.08,500,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',4.5,0.08,1200,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',4.5,0.01,75,45,10)
%
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.02,0.08,250,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.02,0.08,500,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.02,0.08,1200,45,1)
% LVxLMSDeCorrHeterodyne('drwatsonSR8K.wav',0.02,0.01,75,45,10)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

A = real(A);
Mu = real(Mu);
Freq = real(Freq);
NoTaps = floor(abs(NoTaps));
Err = [];

if NoTaps > 200
   NoTaps = 200;
   Comment = 'Limiting Number of LMS Taps to 200'
end

DeCorrDelay = floor(abs(DeCorrDelay));
if DeCorrDelay > 100
   DeCorrDelay = 100;
end

%==========================================================================
% [y,Fs,bits] = wavread('drwatsonSR8K.wav');
[y,Fs,bits] = wavread(strWavFile);
y = y'; lenFile = length(y);
y = (1/max(abs(y)))*y;
t = [0:1:lenFile-1]/Fs;
LMSDeCorDataVec = y + A*sin(2*pi*t*Freq);
SR = 2^ceil(log2(length(LMSDeCorDataVec)))
Fs = Fs
%
ScaleFac = max(abs(LMSDeCorDataVec));
LMSDeCorDataVec = LMSDeCorDataVec/ScaleFac;
TapWt(1:NoTaps,1) = 0;
FiltOut(1,1:lenFile) = 0;
Err(1,1:lenFile) = 0;
n = 1:1:NoTaps;
 
for CurDtaPtr = NoTaps+1:1:lenFile-DeCorrDelay
FiltOut(1,CurDtaPtr) = LMSDeCorDataVec(CurDtaPtr-n)*TapWt;
Err(1,CurDtaPtr)= LMSDeCorDataVec(1,CurDtaPtr + DeCorrDelay-1)-FiltOut(1,CurDtaPtr); 
LIMAdjust = (std(LMSDeCorDataVec(1,CurDtaPtr-(NoTaps):CurDtaPtr)))^2; 
 TapWt = TapWt + 2*Mu*Err(1,CurDtaPtr)/(50*LIMAdjust)*LMSDeCorDataVec(1,CurDtaPtr-n)';
end

figure(17)
clf
% set(17,'color',[1,1,1])

subplot(2,2,1)
plot(LMSDeCorDataVec,'k')
zzc = max(abs(LMSDeCorDataVec));
ylabel(['Amplitude'])
xlabel(['(a) Sample, Input Sig'])
axis([0 inf -zzc zzc])

Err = 0.99*Err/max(abs(Err));

subplot(2,2,3)
plot(Err,'k')
ylabel(['Amplitude'])
xlabel(['(c) Sample, Output/Err'])
axis([0 inf -inf inf])

subplot(2,2,2)
xfft = abs(fft(LMSDeCorDataVec,SR))+(10^-6);
xfft = xfft(1,1:SR/2);
zx = max(xfft);
plot((Fs/2)*[0:1:SR/2-1]/(SR/2),20*log10(xfft/zx),'k')
axis([0, (Fs/2), -inf, inf])
ylabel(['Mag, dB'])
xlabel(['(b) Input Signal Frequency, Hz'])

subplot(2,2,4)
st = fix(length(Err)/2);
theans = abs(fft(Err(1,st:length(Err)),SR))+eps;
theans = theans(1,1:SR/2);
plotlim = max(abs(theans));
plot((Fs/2)*[0:1:SR/2-1]/(SR/2),20*log10(theans/plotlim),'k')
axis([0,(Fs/2),-inf,inf])
ylabel(['Mag, dB'])
xlabel(['(d) Output/Error Frequency, Hz'])

%==========================================================================
ftSR = 4096;
figure(18)
clf
set(18,'color',[1,1,1])
frtw = abs(fft(TapWt,ftSR))';
frtw = frtw/max(frtw);
szfrtw = size(frtw);
hold on
xplot = (Fs/2)*[0:1:ftSR/2]/(ftSR/2);
yplot = (frtw(1,1:ftSR/2+1));
plot(xplot,yplot,'k')
ylabel(['Mag'])
xlabel(['Adaptive Filter Response, Hz'])
axis([0,inf,-inf,inf])
%==========================================================================
%  Linear Plots
figure(177)
clf
% set(17,'color',[1,1,1])

subplot(2,2,1)
plot(LMSDeCorDataVec,'k')
zzc = max(abs(LMSDeCorDataVec));
ylabel(['Amplitude'])
xlabel(['(a) Sample, Input Sig'])
axis([0 inf -zzc zzc])

Err = 0.99*Err/max(abs(Err));

subplot(2,2,3)
plot(Err,'k')
ylabel(['Amplitude'])
xlabel(['(c) Sample, Output/Err'])
axis([0 inf -inf inf])

subplot(2,2,2)
xfft = abs(fft(LMSDeCorDataVec,SR))+(10^-6);
xfft = xfft(1,1:SR/2);
zx = max(xfft);
plot((Fs/2)*[0:1:SR/2-1]/(SR/2),(xfft/zx),'k')
axis([0, (Fs/2), -inf, inf])
ylabel(['Mag'])
xlabel(['(b) Input Signal Frequency, Hz'])

subplot(2,2,4)
st = fix(length(Err)/2);
theans = abs(fft(Err(1,st:length(Err)),SR))+eps;
theans = theans(1,1:SR/2);
plotlim = max(abs(theans));
plot((Fs/2)*[0:1:SR/2-1]/(SR/2),(theans/plotlim),'k')
axis([0,(Fs/2),-inf,inf])
ylabel(['Mag'])
xlabel(['(d) Output/Error Frequency, Hz'])

%--------------------------------------------------------------------------

LMSDeCorDataVec = 0.99*LMSDeCorDataVec/max(abs(LMSDeCorDataVec));

sound(LMSDeCorDataVec,Fs)
pause(6)
sound(Err,8000)
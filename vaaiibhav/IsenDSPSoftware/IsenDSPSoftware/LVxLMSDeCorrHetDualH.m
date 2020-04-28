function LVxLMSDeCorrHetDualH(strWavFile,A,Mu,Freq,NoTaps,DeCorrDelay)
% function LVxLMSDeCorrHetDualH(strWavFile,A,Mu,Freq,NoTaps,DeCorrDelay)
% Has the same input arguments and function as LVxLMSDeCorrHeterodyne, but
% uses a Dual-H architectures for better performance.
% strWavFile is the audio file to open and use as the test signal
% A is the amplitude of an interfering tone (sinusoid) to be added to
% the audio file whose file name or path name (including file name and
% extension) is specified as strWavFile.
% Mu is the LMS update term weight;
% NoTaps is the number of taps to use in the LMS adaptive filter;
% DeCorrDelay is the number of samples to delay the filter input relative to the Plant or channel
% delay.
% LVxLMSDeCorrHeterodyne, but uses a Dual-H architecture.
% The following plots are created: Frequency response (DTFT) of the adaptive filter when converged 
% (magnitude v. frequency in Hz), the test signal (amplitude v. sample),  the filtered test signal (amplitude 
% v. sample), the DTFT of the test signal (magnitude v. frequency in Hz, and the DTFT (magnitude v. 
% frequency in Hz) of the last 50 % of the filtered test signal.
% Test calls:
% LVxLMSDeCorrHetDualH('drwatsonSR8K.wav',0.2,0.03,250,65,3)
% LVxLMSDeCorrHetDualH('drwat8Kplus400HzAnd740Hz.wav',0,0.025,0,73,3)
% LVxLMSDeCorrHetDualH('drwat8Kplus400HzAnd740Hz.wav',0.04,0.025,1150,73,3)
% LVxLMSDeCorrHetDualH('drwatson8KplusFSK400n1100.wav',0,0.025,0,73,3)
% LVxLMSDeCorrHetDualH('drwatson8KplusFSK1305n1465.wav',0,0.025,0,73,3)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

A = real(A);
Mu = real(Mu);
Freq = real(Freq);
NoTaps = floor(abs(NoTaps));
Err = [];

if NoTaps > 400
   NoTaps = 400;
   Comment = 'Limiting Number of LMS Taps to 400'
end

DeCorrDelay = floor(abs(DeCorrDelay));
if DeCorrDelay > 100
   DeCorrDelay = 100;
end

%==========================================================================
% [y,Fs,bits] = wavread('drwatsonSR8K.wav');
[y,Fs,bits] = wavread(strWavFile);
y = y'; 
lenFile = length(y);
y = (1/max(abs(y)))*y;
t = [0:1:lenFile-1]/Fs;
LMSDeCorDataVec = y + A*sin(2*pi*t*Freq);
szLMSDeCorDataVec = size(LMSDeCorDataVec)
SR = 2^ceil(log2(length(LMSDeCorDataVec)))

ScaleFac = max(abs(LMSDeCorDataVec));
LMSDeCorDataVec = LMSDeCorDataVec/ScaleFac;
TapWt(1:NoTaps,1) = 0;
FiltOut(1,1:lenFile) = 0;
Err(1,1:lenFile) = 0;
n = 1:1:NoTaps;
%=========================================================================
BestTapWt(1,1:NoTaps) = 0;
TestTapWt = zeros(1,NoTaps);
ActualErr = zeros(1,lenFile);
TestErr = zeros(1,lenFile);
TestErr(1,1:NoTaps) = 1;
CurBestERLE(1,lenFile) = zeros;
CurBestERLE(1,1:NoTaps) = 1;
TestERLE = 0.001;
Start = NoTaps + 1;

%=========================================================================
 
for CurDtaPtr = NoTaps+1:1:lenFile-DeCorrDelay
    
if TestERLE > CurBestERLE(1,CurDtaPtr-1)
   BestTapWt = TestTapWt;
   CurBestERLE(1,CurDtaPtr) = TestERLE;
else
   CurBestERLE(1,CurDtaPtr) = CurBestERLE(1,CurDtaPtr-1);
end             
ActualFiltOut =  sum(BestTapWt.*LMSDeCorDataVec(1,CurDtaPtr-n));
TestFiltOut =  sum(TestTapWt.*LMSDeCorDataVec(1,CurDtaPtr-n));

lmsd = LMSDeCorDataVec(1,CurDtaPtr + DeCorrDelay-1);
Err(1,CurDtaPtr) = lmsd - ActualFiltOut;          
TestErr(1,CurDtaPtr) = lmsd - TestFiltOut;  
TestERLE = sum(LMSDeCorDataVec(1,CurDtaPtr-n).^2)/(sum(TestErr(1,CurDtaPtr-n).^2)+10^-8);  
%TestERLE = sum(LMSDeCorDataVec(1,CurDtaPtr).^2)/(sum(TestErr(1,CurDtaPtr).^2)+10^-8);  
TestTapWt = TestTapWt + ...
   2*Mu*TestErr(1,CurDtaPtr)*(LMSDeCorDataVec(1,CurDtaPtr-n))/sum((LMSDeCorDataVec(1,CurDtaPtr-n)).^2 + 10^-8);
end

figure(171)
clf

subplot(2,1,1)
plot(LMSDeCorDataVec)
zzc = max(abs(LMSDeCorDataVec));
ylabel(['Amplitude'])
xlabel(['(a) Sample, Input Signal'])
axis([0 inf -zzc zzc])

Err = 0.99*Err/max(abs(Err));

subplot(2,1,2)
plot(Err)
vcx = 1;
ylabel(['Amplitude'])
xlabel(['(b) Sample, Output/Error'])
axis([0 inf -inf inf])

figure(172)
clf
subplot(2,1,1)
szLMSDeCorDataVec = size(LMSDeCorDataVec)
SR = SR
xfft = abs(fft(LMSDeCorDataVec,SR));
zx = max(xfft);
plot((Fs/2)*[0:1:SR/2]/(SR/2),( xfft(1,1:(SR/2 + 1))/zx)  )
axis([0, Fs/2, -inf, inf])
ylabel(['Mag'])
xlabel(['(a) Frequency, Hz, Test Signal'])

subplot(2,1,2)
st = fix(length(Err)/2);
theans = abs(fft(Err(1,st:length(Err)),SR));
plotlim = max(abs(theans));
plot((Fs/2)*[0:1:SR/2]/(SR/2),(theans(1,1:(SR/2 + 1))/plotlim))
axis([0,(Fs/2),-inf,inf])
ylabel(['Mag'])
xlabel(['(b) Frequency, Hz, Filtered Output Signal'])

ftSR = 4096;
figure(181)
clf
subplot(211)
stem(BestTapWt)
xlabel('(a) Adaptive Filter Tap Number')
ylabel('Weight')

subplot(212)
frtw = abs(fft(BestTapWt,ftSR));
frtw = frtw/max(frtw);
szfrtw = size(frtw);
hold on
xplot = (Fs/2)*[0:1:ftSR/2]/(ftSR/2);
yplot = (frtw(1,1:ftSR/2+1));
plot(xplot,yplot)
ylabel(['Mag'])
xlabel(['(b) Adaptive Filter Response, Hz'])
axis([0,inf,-inf,inf])

figure(221)
clf
plot(20*log10(CurBestERLE(1,1:length(CurBestERLE))))
plotmin = 1.2*min(20*log10(CurBestERLE(1,1:length(CurBestERLE))));
plotmax = 1.2*max(20*log10(CurBestERLE(1,1:length(CurBestERLE))));
ylabel(['Mag, dB'])
xlabel('Sample')

sound(LMSDeCorDataVec,Fs)
pause(6)
sound(Err,8000)


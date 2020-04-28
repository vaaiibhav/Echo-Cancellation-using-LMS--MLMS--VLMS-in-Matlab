function LVxInterp8Kto11025(testSigType,tsFreq)
% LVxInterp8Kto11025(testSigType,tsFreq)
% Converts an audio signal having a sample rate of 8 kHz to one sampled at
% 11.025 kHz using both sinc and linear interpolation. The audio files used
% can be either the audio file 'drwatsonSR8K.wav' or a cosine of user 
% designated frequency. For sinc interpolation it uses 10 samples of the input 
% signal and a 100 by 10 sinc interpolation matrix to generate sample values 
% located at the fractional sample index values 1:320/441:length(OriginalAudioFile)
% To use 'drwatsonSR8K.wav', pass testSigType as 1 and tsFreq as [] or any
% number; to use a cosine, pass testSigType as 0 and pass tsFreq as the
% desired cosine frequency in Hz.
% Test calls:
% LVxInterp8Kto11025(0,250)
% LVxInterp8Kto11025(0,1950)
% LVxInterp8Kto11025(1,[])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
N = 1000;
SampDecRate = 100;
n = -N/2:1:N/2;
n = n/SampDecRate;
ysinc = sinc(n); 
ysinc = (hamming(length(ysinc))').*ysinc;
sincMat = [];
%==========================================================================
sincMat([2*N - SampDecRate],fix(N/SampDecRate)) = 0; 
for ctr = 1:1:fix(N/SampDecRate)
sincMat((ctr-1)*SampDecRate + 1:(ctr-1)*SampDecRate + N + 1,ctr) = ysinc';
end
%==========================================================================
netSincMat = sincMat(N-SampDecRate:N,:);
sznetSincMat = size(netSincMat);
NoSamps = sznetSincMat(2);
sincMat = [];

if testSigType==1
strAudio = 'drwatsonSR8K.wav';
[testSig,FS,NBITS]= wavread(strAudio);
testSig = testSig';
lenTestSig = length(testSig);
durtestSig = lenTestSig/FS;
else
  t = 0:1/7999:1;
  testSig = cos(2*pi*t*tsFreq);
  FS = 8000;
  lenTestSig = 8000;
  durtestSig = 1;
end
    
RatNum = 320;
RatDen = 441;

newIndices = ([1:RatNum/RatDen:lenTestSig]');

indTable = zeros(length(newIndices),5);

indTable(1:length(newIndices),1) = newIndices;
indTable(1:length(newIndices),2) = floor(newIndices);
indTable(1:length(newIndices),3) = round(100*(newIndices-floor(newIndices)))+1;
szindTable = size(indTable);

LoLim = ceil(NoSamps/2)+ 2;
if testSigType==0
UpLim = 2100; 
elseif testSigType==1
UpLim = szindTable(1)-(ceil(NoSamps/2)+ 2);
end
InterpSamp = zeros(UpLim,1);
LinInterpSamp = zeros(UpLim,1);

for Ctr = LoLim:1:UpLim
%    theCtr = Ctr
    xx = indTable(Ctr,2:3);
    RowNo = xx(2);
    SampLo = xx(1);
    yy = testSig(SampLo-4:SampLo+5);
indTable(Ctr,4) = netSincMat(RowNo,:)*yy'; %sinc interp
   LinDel = (RowNo-1)/100;
indTable(Ctr,5) = yy(5)*(1-LinDel) + LinDel*yy(6); % lin interp
end

[b,a] = cheby1(12,0.5,7.9/11.025);
SincInterpSig = filter(b,a,indTable(:,4));
LinInterpSig = filter(b,a,indTable(:,5));

figure(88)
lenIntSamp = length(newIndices);
plotlim = 1.2*max(abs(testSig(1,1:1105)));
subplot(311) 
plot(testSig) 
xlabel('(a) Sample, Test Signal (SR = 8 kHz)')
ylabel('Amp')
axis([0,800,-plotlim,plotlim])

subplot(312)
plot(indTable(:,4))
xlabel('(b) Sample, Sinc-Interpolation (SR = 11.025 kHz)')
ylabel('Amp')
axis([0,1103,-plotlim,plotlim])

subplot(313)
plot(indTable(:,5))
xlabel('(c) Sample, Linear Interpolation (SR = 11.025 kHz)')
ylabel('Amp')
axis([0,1103,-plotlim,plotlim])

n = 1500;
figure(86)
clf
hold on
LoLim = n+1;

HiLim = n+10;
xvec = indTable(LoLim:HiLim,1);

disTestSig = indTable(LoLim:HiLim,2);
minTSigind = min(disTestSig);
maxTestSigind = max(disTestSig);

stem([minTSigind:maxTestSigind],testSig(minTSigind:maxTestSigind),'bo')
stem(xvec,indTable(LoLim:HiLim,4),'k*')
stem(xvec,indTable(LoLim:HiLim,5),'rd')
xlabel('Input Sequence Sample')
ylabel('Amplitude')
axis([minTSigind-0.5,maxTestSigind+0.5,-inf,inf])

indvec = 1001:1001+1023;
lenindvec = length(indvec)
% plot raw signal and interpolated spectra (i.e. prior to post inter lp
% filtering)
s0 = testSig(indvec)'.*hamming(1024);
s1 = indTable(indvec,4).*hamming(1024);
s2 = indTable(indvec,5).*hamming(1024);

f0 = abs(fft(s0,2^17));
f0 = f0/max(f0);
f0 = 20*log10(f0+eps);

f1 = abs(fft(s1,2^17));
f1 = f1/max(f1);
f1 = 20*log10(f1+eps);
f2 = abs(fft(s2,2^17));
f2 = f2/max(f2);
f2 = 20*log10(f2+eps);

figure(998)
clf

subplot(131)
plot([0:2^16-1]/(2^16),f0(1:2^16))
xlabel('(a) Norm Freq')
ylabel('Mag, dB')
grid on
axis([0,1,-120,10])

subplot(132)
plot([0:2^16-1]/(2^16),f1(1:2^16))
xlabel('(b) Norm Freq')
grid on
ylabel('Mag, dB')
axis([0,1,-120,10])

subplot(133)
plot([0:2^16-1]/(2^16),f2(1:2^16))
xlabel('(c) Norm Freq')
grid on
ylabel('Mag, dB')
axis([0,1,-120,10])

% now compute & display spectra after post-interpolation lp filtering
% plot raw signal and interpolated spectra (i.e. prior to post-interp
% lowpass filtering)

s1 = SincInterpSig(indvec).*hamming(1024);
s2 = LinInterpSig(indvec).*hamming(1024);

f1 = abs(fft(s1,2^17));
f1 = f1/max(f1);
f1 = 20*log10(f1+eps);
f2 = abs(fft(s2,2^17));
f2 = f2/max(f2);
f2 = 20*log10(f2+eps);

figure(1998)
clf

subplot(131)
plot([0:2^16-1]/(2^16),f0(1:2^16))
xlabel('(a) Norm Freq')
ylabel('Mag, dB')
grid on
axis([0,1,-120,10])

subplot(132)
plot([0:2^16-1]/(2^16),f1(1:2^16))
xlabel('(b) Norm Freq')
grid on
ylabel('Mag, dB')
axis([0,1,-120,10])

subplot(133)
plot([0:2^16-1]/(2^16),f2(1:2^16))
xlabel('(c) Norm Freq')
grid on
ylabel('Mag, dB')
axis([0,1,-120,10])


sound(testSig,FS)  
pause(1.2*durtestSig)
sound(indTable(:,4),11025)
pause(1.2*durtestSig)
sound(indTable(:,5),11025)

pause(2)
sound(testSig,FS)  
pause(1.2*durtestSig)
sound(SincInterpSig,11025)
pause(1.2*durtestSig)
sound(LinInterpSig,11025)




    
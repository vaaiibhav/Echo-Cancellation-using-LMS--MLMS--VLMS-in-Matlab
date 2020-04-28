function [ToneFreq,Fnyq,BnSp] = LVx_DetectContTone(A,Freq,RorSS,SzWin,OvrLap,AudSig)
% function [ToneFreq,Fnyq,BnSp] = LVx_DetectContTone(A,Freq,RorSS,SzWin,OvrLap,AudSig)
% Mixes a tone of amplitude A and frequency Freq with an audio file and
% attempts to identify the frequency of the interfering tone, which may
% have a steady-state amplitude A when RorSS is passed as 0, or has an
% amplitude that linearly ramps from 0 to A over the length of the audio
% file when RorSS is passed as 1. The audio file is 'drwatsonSR8K.wav', 
% 'whoknowsSR8k.wav', or white noise, which are selected respectively by passing
% AudSig as 1, 2, or 3. SzWin is the size of time window in
% samples into which the test signal (audio file plus tone) is partitioned
% for analysis. In partitioning the test signal into time windows or frames, an
% overlap of either zero samples or fix(SzWin/2) is performed, according to whether
% OvrLap is passed as 0 or 1, respectively.
% If OvrLap is paased as 0, a hamming window is applied to the frames for
% frequency analysis, but not directly to the frames of TDMat. If OvrLap is passed
% as 1, meaning a 50% sample overlap is called for, the frames of TDMat are
% windowed with a hamming window, and the frequency analysis is performed
% by taking the DFT of the frames of TDMat.
% A number of figures are created,
% including a first figure displaying mean bin magnitude of all bins for
% each frame, which serves to identify periods of high and low energy in
% the signal, corresponding to active speech and background sound, second
% and third figures that are 3-D spectrogram plots of DFT magnitude versus Bin and
% Frame for the complete signal and a version based on only the background or
% relatively low energy frames. Fourth and fifth figures display the
% normalized bin derivatives versus frame for the complete signal and the low
% portion of the signal.
% The output arguments consist of a list of possible interfering 
% (relatively steady-state) tones in Hz as ToneFreq, the Nyquist rate as Fnyq, and
% the bin spacing in Hz. This allows another program to determine if a list of
% frequencies provided as ToneFreq contains adjacent bin frequencies and are 
% therefor likely to be members of the same spectral component. 
% Test calls
%
% [T,F,B] = LVx_DetectContTone(0.011,100,0,512,1,1)
% [T,F,B] = LVx_DetectContTone(0.008,100,0,512,0,1)
% [T,F,B] = LVx_DetectContTone(0.005,94,0,512,0,1)
% [T,F,B] = LVx_DetectContTone(0.07,150,0,512,1,3)
%
% [T,F,B] = LVx_DetectContTone(0.011,200,0,512,1,1)
% [T,F,B] = LVx_DetectContTone(0.009,200,0,512,0,1)
% [T,F,B] = LVx_DetectContTone(0.008,203.125,0,512,0,1)
% [T,F,B] = LVx_DetectContTone(0.01,200,0,512,1,2)
% [T,F,B] = LVx_DetectContTone(0.1,200,0,512,1,3)
% [T,F,B] = LVx_DetectContTone(0.01,210,0,512,1,1)
% [T,F,B] = LVx_DetectContTone(0.025,200,1,512,1,1)
% [T,F,B] = LVx_DetectContTone(0.02,200,1,512,1,2)
% [T,F,B] = LVx_DetectContTone(0.1,200,1,512,1,3)
%
% [T,F,B] = LVx_DetectContTone(0.005,500,0,512,1,1)
% [T,F,B] = LVx_DetectContTone(0.007,500,0,512,0,1)
% [T,F,B] = LVx_DetectContTone(0.01,500,0,512,1,2)
% [T,F,B] = LVx_DetectContTone(0.08,500,0,512,1,3)
% [T,F,B] = LVx_DetectContTone(0.02,500,1,512,1,1)
% [T,F,B] = LVx_DetectContTone(0.018,500,1,512,1,2)
% [T,F,B] = LVx_DetectContTone(0.1,500,1,512,1,3)
% [T,F,B] = LVx_DetectContTone(0.005,1000,0,512,0,1)
% [T,F,B] = LVx_DetectContTone(0.004,1000,0,512,1,1)
%
% [T,F,B] = LVx_DetectContTone(0.002,3000,0,512,1,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
fr = [];

if AudSig==1
[y,Fs,bits] = wavread('drwatsonSR8K.wav');
elseif AudSig==2
[y,Fs,bits] = wavread('whoknowsSR8K.wav');
elseif AudSig==3
Fs = 8000;
y = randn(5*Fs,1);
else
  [y,Fs,bits] = wavread('drwatsonSR8K.wav');  
end

NyqRate = Fs/2;
Fnyq = NyqRate;
y = y'; 
lenFile = length(y);
y = (1/max(abs(y)))*y;
y = y - mean(y);
t = [0:1:lenFile-1]/Fs;
if RorSS==0
  TstSig = y + A*cos(2*pi*t*Freq);
else
  TstSig = y + ([1:1:lenFile]/lenFile).*(A*cos(2*pi*t*Freq));
end
if ~(rem(SzWin,2)==0)
    error('SzWin must be an even number! Exiting...')
    return
end

n = 2^(ceil(log2(lenFile)))
Fr = abs(fft(TstSig,n));
figure(188)
xvec = (Fs/2)*[0:1:length(Fr)/2]/(length(Fr)/2);
plot(xvec,Fr(1,1:fix(length(Fr)/2)+1))
xlabel('Frequency, Hz')
ylabel('Magnitude')

if OvrLap==0
TDMat = LVxVector2FramesInMatrix(TstSig,SzWin,0);
FDMat4Recon = fft(TDMat);
szy = size(TDMat);
lgFty = fft((hamming(SzWin)*ones(1,szy(2))).*TDMat);
SampOverLap = 0;
else
SampOverLap = fix(SzWin/2);
TDMat = LVxVector2FramesInMatrix(TstSig,SzWin,SampOverLap);
szy = size(TDMat);
TDMat = (hamming(SzWin)*ones(1,szy(2))).*TDMat;
lgFty = fft(TDMat);
end

fty = lgFty(1:SzWin/2+1,:);
szfty = size(fty);

afty = abs(fty);
afty = afty/(max(max(afty)));

%==========================================================================
meanframes = mean(afty);
[Nhist,xhist] = hist(meanframes);
percentlow = Nhist(1)/sum(Nhist);
minFrameMat = [];
histthresh = 1;
histthresh = histthresh-1;
framecnt = 1;
while framecnt < 15
histthresh = histthresh + 1;
[i,j] = find(meanframes < xhist(histthresh));

minFrameMat = afty(:,j);
szminFrameMat = size(minFrameMat);

meanmag = mean(meanframes);
logfabsfty = 20*log10(afty+10^(-8));
logmeanmag = 20*log10(meanmag);

meanMinMat = mean(minFrameMat);
meanMinMag = mean(meanMinMat);
logMinMat = 20*log10(minFrameMat+10^(-8));
logmeanMinMag = 20*log10(meanMinMag);
szlogminmat = size(logMinMat);
framecnt = szlogminmat(2);
end

figure(122)
clf
hold on
stem(meanframes)
hline = line([0,length(meanframes)],[xhist(histthresh),xhist(histthresh)]);
set(hline,'linestyle',':')
xlabel('Frame (Entire Signal)')
ylabel('Avg Bin Mag')

figure(98)
clf
szlogabsfty = size(logfabsfty);
xvec = NyqRate*[0:1:szlogabsfty(1)-1]/szlogabsfty(1);
yvec = [0:1:szlogabsfty(2)-1];
mesh(yvec,xvec,logfabsfty);
xlabel('Frame (Entire Signal)')
ylabel('Frequency, Hz')
zlabel('Mag, dB')

figure(99)
clf

xvec = NyqRate*[0:1:szlogminmat(1)-1]/szlogminmat(1);
yvec = [0:1:szlogminmat(2)-1];

mesh(yvec,xvec,logMinMat)
xlabel('Frame (Min Mag Frames)')
ylabel('Frequency, Hz')
zlabel('Mag, dB')

evalMat = zeros(szfty);
hiAmps = find(logfabsfty >logmeanmag+3);
evalMat(hiAmps) = 1;
sumrows = sum(evalMat');
totalFrames = szfty(2); % 
[N,X] = hist(sumrows);
hisumrows = X(10);
hiBINS = find(sumrows>=hisumrows); 
intFREQhiAmpFull = (Fs/2)*(hiBINS-1)/(fix(SzWin/2))

% =========================deriv===========================================
deriv = sqrt((1/(szfty(2)-1))*(afty(:,2:szfty(2)) - afty(:,1:szfty(2)-1)).^2);
meanrows = mean(afty');
avgDeriv = mean(deriv')./meanrows;
[Nder,Xderiv] = hist(avgDeriv);
loDerivThreshFull = Xderiv(3);
loDerivIndsFull = find(avgDeriv <= loDerivThreshFull);
intFREQDerivFull = (Fs/2)*(loDerivIndsFull-1)/(fix(SzWin/2))
%==========================================================================
sgnderiv = afty(:,2:szfty(2)) - afty(:,1:szfty(2)-1);
posderivF = find(sgnderiv > 0);
zeroderivF = find(sgnderiv == 0);
sgnderiv = -ones(size(sgnderiv));
sgnderiv(posderivF) = 1;
sgnderiv(zeroderivF) = 0;
sumrowsposderiv = sum(sgnderiv');
stdDeriv = std(sumrowsposderiv);
hiPosDerivInds = find(sumrowsposderiv>=3*stdDeriv);
intFREQhiPosDerivAllFrames = (Fs/2)*(hiPosDerivInds-1)/(fix(SzWin/2))
%==========================================================================
figure(554)
clf
xvec = [0:1:length(avgDeriv)-1]/(length(avgDeriv)-1)*(Fs/2);
stem(xvec,avgDeriv)
xlabel('Frequency, Hz')
ylabel('Avg Mag F.O.D. for all Frames')

% do with MinFrameMat
evalMinMat = zeros(szminFrameMat);
hiAmpsMinMat = find(logMinMat >logmeanMinMag+8);
evalMinMat(hiAmpsMinMat) = 1;
sumrowsMM = sum(evalMinMat');
MinMatTotalFrames = szminFrameMat(2);
% pick multiple top candidates
[Nmm,Xmm] = hist(sumrowsMM);
hisumrowsthr = Xmm(10);
hiBINSmm = find(sumrowsMM>=hisumrowsthr); 
interferFREQSmm = (Fs/2)*(hiBINSmm-1)/(fix(SzWin/2))

% deriv MinMat=============================================================
magderivMM = sqrt((1/(szminFrameMat(2)-1))*(minFrameMat(:,2:szminFrameMat(2)) - minFrameMat(:,1:(szminFrameMat(2)-1))).^2);
meanrowsMM = mean(minFrameMat');
normDerivMM = (mean(magderivMM'))./meanrowsMM;
[Ndermm,Xderivmm] = hist(normDerivMM);
loDerivThresh = Xderivmm(3);
loDerivInds = find(normDerivMM <= loDerivThresh);
interFREQsDerivMinMag = (Fs/2)*(loDerivInds-1)/(fix(SzWin/2))

meannormDerivMM = mean(normDerivMM);
stdnormDerivMM = std(normDerivMM);
supLoDerivMM = find(normDerivMM <= (meannormDerivMM-3*stdnormDerivMM));
supLoDerivFreqsMinMag = (Fs/2)*(supLoDerivMM-1)/(fix(SzWin/2))

figure(555)
clf
xvec = ([0:1:length(normDerivMM)-1]/(length(normDerivMM)-1))*(Fs/2);
stem(xvec,normDerivMM)
xlabel('Frequency, Hz')
ylabel('Avg Mag F.O.D. for low Mag Frames')

sgnderivMM = minFrameMat(:,2:szminFrameMat(2)) - minFrameMat(:,1:(szminFrameMat(2)-1));
posderiv = find(sgnderivMM > 0);
zeroderiv = find(sgnderivMM == 0);
sgnderivMM = -ones(size(sgnderivMM));
sgnderivMM(posderiv) = 1;
sgnderivMM(zeroderiv) = 0;
sumrowsposderiv = sum(sgnderivMM');
avgDerivMM = mean(sumrowsposderiv);
stdDerivMM = std(sumrowsposderiv);
hiPosDerivInds = find(sumrowsposderiv>=3*stdDerivMM);
interFREQshiPosDeriv = (Fs/2)*(hiPosDerivInds-1)/(fix(SzWin/2))

%interFreqMinMatHiMag
InterferSinusoids = [];
lenHiMagsVec = length(interferFREQSmm);
for ctr = 1:1:lenHiMagsVec
    y = find(interFREQsDerivMinMag-interferFREQSmm(ctr)==0);
    if ~(isempty(y))
        InterferSinusoids = [InterferSinusoids,interferFREQSmm(ctr)];
    end
end

%interFREQsDerivMinMag
if ~(isempty(interFREQshiPosDeriv))
    for ctr = 1:1:lenHiMagsVec
        y = find(interFREQshiPosDeriv-interferFREQSmm(ctr)==0);
        if ~(isempty(y))
        InterferSinusoids = [InterferSinusoids,interferFREQSmm(ctr)];
        end
    end
    for ctr = 1:1:length(interFREQsDerivMinMag)
        y = find(interFREQshiPosDeriv-interFREQsDerivMinMag(ctr)==0);
        if ~(isempty(y))
        InterferSinusoids = [InterferSinusoids,interFREQsDerivMinMag(ctr)];
        end
    end   
end
BnSp = (Fs/2)/fix(SzWin/2);
ToneFreq = [InterferSinusoids,supLoDerivFreqsMinMag];
ToneFreq = sort(ToneFreq);
% eliminate any duplicates
lenTF = length(ToneFreq);
x = ToneFreq(1,2:lenTF) - ToneFreq(1,1:lenTF-1);
xzzs = find(x==0);
ToneFreq(xzzs) = [];
%=========================

sound(TstSig,8000) %===============================================

 if OvrLap==0
     ProperFDMat = FDMat4Recon;
 else
     ProperFDMat = lgFty;
 end
% reconstruct signal frames==========================================
NewTD = real(ifft(ProperFDMat))';
tdSig = LVxTDMat2SigVec(NewTD,SzWin,SampOverLap);

pause(5)
sound(tdSig,8000)





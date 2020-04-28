function [ToneFreq,Fnyq,BnSp] = LVx_DetnFiltFIRContTone(A,Freq,RorSS,SzWin,OvrLap,AudSig,Rp,As)
% function [ToneFreq,Fnyq,BnSp] = LVx_DetnFiltFIRContTone(A,Freq,RorSS,SzWin,OvrLap,AudSig,Rp,As)
% Creates a test signal comprising an audio signal and an interfering
% sinusoid, and attempts to identify the interfering tone. The input and
% output arguments are identical to those of the script ML_DetectContTone
%
% This script then evaluates the list of candidate interfering tones
% ToneFreq to determine if filtering should be performed. Frequencies below 80
% Hz are eliminated from the list as the two audio files contain prominent
% 60 Hz components, and the program is only designed to filter out one main
% spectral component which lies above 80 Hz. Also, the ear's
% sensitivity at low frequencies is very low, so the audibility of such
% tones is not high, and is often well-masked by other signal components.
% The candidate frequencies above 80 Hz comprise either a single frequency or 
% a number of contiguous frequencies(i.e.,lying in adjacent bins of the
% DFT). When a number of contiguous candidate frequencies exist, the lower
% and upper frequency bounds are established. If only a single frequency
% candidate exists, upper and lower bounds surrounding it are created so
% that a notch filter of reasonable width can be designed. When the
% interfering frequency or frequencies lie at or below 0.035pi radians 
% (about 140 Hz for this exercise), a highpass filter is used; when the 
% interfering frequencies lie above 0.95pi radians, a lowpass filter is used.
% Otherwise, a notch filter is used.
% The equiripple FIR filter is designed to have maximum passband ripple of Rp dB and
% minimum stopband attenation of As dB.
% Filtering is performed two ways, first, by time domain convolution of the
% entire test signal with the designed filter, and secondly, by frequency domain convolution 
% that uses the original analysis matrix, which is doubled in column length 
% (by addition of zeros to each column) to accomodate both the analysis
% and frequency domain convolution. After the columns of this special
% matrix ahve been multiplied by the DFT of the designed filter's impulse response, 
% the IDFT is obtained to form filtered time domain signal frames, and the 
% output signal vector is reconstructed by concatenating frames with the
% proper overlap, if any.
% The test signal, the time domain filtered version, and the frequency
% domain filtered version are all played out through the computer's audio
% system for comparison.
% Several figures are created in addition to those associated with frequency detection per se, 
% including the filter impulse response, the filter frequency response, and the test signal frequency
% content after being filtered as one signal vector, and after being filtered using the DFT Convolution 
% method, reconstructing the signal vector from frames 
%
% [T,F,B] = LVx_DetnFiltFIRContTone(0.015,200,0,512,1,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.01,210,0,512,1,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.01,200,0,512,1,2,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.1,200,0,512,1,3,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.01,500,0,512,1,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.01,500,0,512,1,2,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.08,500,0,512,1,3,1,30)
%
% [T,F,B] = LVx_DetnFiltFIRContTone(0.025,200,1,512,1,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.02,200,1,512,1,2,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.1,200,1,512,1,3,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.02,500,1,512,1,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.018,500,1,512,1,2,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.1,500,1,512,1,3,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.07,150,0,512,1,3,1,30)
%
% [T,F,B] = LVx_DetnFiltFIRContTone(0.004,1000.5,0,512,0,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.005,1000.5,0,512,0,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.004,1000.5,0,512,1,1,1,30)
% [T,F,B] = LVx_DetnFiltFIRContTone(0.005,1000.5,0,512,1,1,1,30)
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
TstSig = TstSig*(1/max(abs(TstSig)));

n = 2^(ceil(log2(lenFile)))
Fr = abs(fft(TstSig,n));

figure(188)

xvec = (Fs/2)*[0:1:length(Fr)/2]/(length(Fr)/2);
plot(xvec,Fr(1,1:fix(length(Fr)/2)+1))
xlabel('Frequency, Hz (Test Signal)')
ylabel('Magnitude')

if OvrLap==0
TDMat = LVxVector2FramesInMatrix(TstSig,SzWin,0);
szy = size(TDMat);
WinTDMat = (hamming(szy(1))*ones(1,szy(2))).*TDMat;
dblWinTDMat = WinTDMat; % dbls length of cols, 2nd half of cols = zeros
dblWinTDMat(szy(1)+1:2*szy(1),:) = 0;
lgFty = fft(dblWinTDMat);
dblTDMat = TDMat; % dbls length of cols, 2nd half of cols = zeros
dblTDMat(szy(1)+1:2*szy(1),:) = 0;
FDMat4Recon = fft(dblTDMat);
SampOverLap = 0;
else
SampOverLap = fix(SzWin/2);
TDMat = LVxVector2FramesInMatrix(TstSig,SzWin,SampOverLap);
szy = size(TDMat);
TDMat = hamming(SzWin)*ones(1,szy(2)).*TDMat;
dblTDMat = TDMat; % doubles length of columns, second half of cols = zeros
dblTDMat(szy(1)+1:2*szy(1),:) = 0;
Sszy = size(dblTDMat);
lgFty = fft(dblTDMat,2*SzWin);
end

%lgFty = fft(dblTDMat,2*SzWin);
szlgFty = size(lgFty)

fty = lgFty(1:(szlgFty(1)/2+1),:); % pos bins only for frequency analysis
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


stem(meanframes)
hold on
hline = line([0,length(meanframes)],[xhist(histthresh),xhist(histthresh)]);
hold off
%set(hline,'linestyle',':')
xlabel('Frame (Entire Signal)')
ylabel('Avg Bin Mag')

figure(98)

mesh(logfabsfty)
xlabel('Frame (Entire Signal)')
ylabel('Bin')
zlabel('Mag, dB')

figure(99)

mesh(logMinMat)
xlabel('Frame (Min Mag Frames)')
ylabel('Bin')
zlabel('Mag, dB')

evalMat = zeros(szfty);
hiAmps = find(logfabsfty >logmeanmag+3);
evalMat(hiAmps) = 1;
sumrows = sum(evalMat');
totalFrames = szfty(2); % 

[N,X] = hist(sumrows);
hisumrows = X(10);
hiBINS = find(sumrows>=hisumrows); 
intFREQhiAmpFull = (Fs/2)*(hiBINS-1)/(fix(SzWin))

% =========================deriv===========================================
deriv = sqrt((1/(szfty(2)-1))*(afty(:,2:szfty(2)) - afty(:,1:szfty(2)-1)).^2);
meanrows = mean(afty');
avgDeriv = mean(deriv')./meanrows;
[Nder,Xderiv] = hist(avgDeriv);
loDerivThreshFull = Xderiv(3);
loDerivIndsFull = find(avgDeriv <= loDerivThreshFull);
intFREQDerivFull = (Fs/2)*(loDerivIndsFull-1)/(fix(SzWin))
%==========================================================================
sgnderiv = afty(:,2:szfty(2)) - afty(:,1:szfty(2)-1);
posderivF = find(sgnderiv > 0);
zeroderivF = find(sgnderiv == 0);
sgnderiv = -ones(size(sgnderiv));
sgnderiv(posderivF) = 1;
sgnderiv(zeroderivF) = 0;
sumrowsposderiv = sum(sgnderiv');
%avgDerivF = mean(sumrowsposderiv);
stdDeriv = std(sumrowsposderiv);
hiPosDerivInds = find(sumrowsposderiv>=3*stdDeriv);
intFREQhiPosDerivAllFrames = (Fs/2)*(hiPosDerivInds-1)/(fix(SzWin))
%==========================================================================
figure(554)
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
interferFREQSmm = (Fs/2)*(hiBINSmm-1)/(fix(SzWin))

% deriv MinMat=============================================================
magderivMM = sqrt((1/(szminFrameMat(2)-1))*(minFrameMat(:,2:szminFrameMat(2)) - minFrameMat(:,1:szminFrameMat(2)-1)).^2);
meanrowsMM = mean(minFrameMat');
normDerivMM = mean(magderivMM')./meanrowsMM;
[Ndermm,Xderivmm] = hist(normDerivMM);
loDerivThresh = Xderivmm(3);
loDerivInds = find(normDerivMM <= loDerivThresh);
interFREQsDerivMinMag = (Fs/2)*(loDerivInds-1)/(fix(SzWin))

meannormDerivMM = mean(normDerivMM);
stdnormDerivMM = std(normDerivMM);
supLoDerivMM = find(normDerivMM <= (meannormDerivMM-3*stdnormDerivMM));
supLoDerivFreqsMinMag = (Fs/2)*(supLoDerivMM-1)/(fix(SzWin))

figure(555)
xvec = [0:1:length(normDerivMM)-1]/(length(normDerivMM)-1)*(Fs/2);
stem(xvec,normDerivMM)
xlabel('Frequency, Hz')
ylabel('Avg Mag F.O.D. for low Mag Frames')

sgnderivMM = minFrameMat(:,2:szminFrameMat(2)) - minFrameMat(:,1:szminFrameMat(2)-1);
posderiv = find(sgnderivMM > 0);
zeroderiv = find(sgnderivMM == 0);
sgnderivMM = -ones(size(sgnderivMM));
sgnderivMM(posderiv) = 1;
sgnderivMM(zeroderiv) = 0;
sumrowsposderiv = sum(sgnderivMM');
avgDerivMM = mean(sumrowsposderiv);
stdDerivMM = std(sumrowsposderiv);
hiPosDerivInds = find(sumrowsposderiv>=3*stdDerivMM);
interFREQshiPosDeriv = (Fs/2)*(hiPosDerivInds-1)/(fix(SzWin))

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
BnSp = (Fs/2)/SzWin;
ToneFreq = [InterferSinusoids,supLoDerivFreqsMinMag];
ToneFreq = sort(ToneFreq);
% eliminate any duplicates
lenTF = length(ToneFreq);
x = ToneFreq(1,2:lenTF) - ToneFreq(1,1:lenTF-1);
xzzs = find(x==0);
ToneFreq(xzzs) = [];
%=========================

NToneFreq = ToneFreq/Fs; % normalized interfering frequencies

sound(TstSig,8000)

% Filter out a selected frequency or group of frequencies.
FiltCands = ToneFreq;
loFrqs = find(FiltCands < 80); % don't consider low frequencies
FiltCands(loFrqs) = [];
lenFiltCands = length(FiltCands)
if lenFiltCands > 1
frDiffs = FiltCands(1,2:lenFiltCands) - FiltCands(1,1:lenFiltCands-1);
adjFrqs = find(frDiffs==BnSp);
else
    adjFrqs = [];
end

if ~(isempty(adjFrqs)) 
lenadjFrqs = length(adjFrqs);
adjFrqs = [adjFrqs,(adjFrqs(lenadjFrqs)+1)];
lenadjFrqs = length(adjFrqs);
maxFrqDiff = (FiltCands(adjFrqs(lenadjFrqs)) - FiltCands(adjFrqs(1)))
    if maxFrqDiff <= 50
        DoFilt = 1;
        comment = 'max freq diff <= 50'
    else
        DoFilt = 0;
        comment = 'max freq diff > 50'
    end
else
 if length(FiltCands)==0
    DoFilt=0;
 else 
 FiltCands = FiltCands(length(FiltCands));
 DoFilt = 1;
 end
end
  
if DoFilt==1

 if length(FiltCands)>1
    minNormNotchFreq = FiltCands(adjFrqs(1))/Fnyq
    maxNormNotchFreq = FiltCands(adjFrqs(lenadjFrqs))/Fnyq
    if maxNormNotchFreq - minNormNotchFreq < 0.015;
        MMmean = (maxNormNotchFreq + minNormNotchFreq)/2;
        minNormNotchFreq = MMmean - 0.0075;
        maxNormNotchFreq = MMmean + 0.0075;
    end
 else
     minNormNotchFreq = FiltCands(1)/Fnyq - 0.0075
     maxNormNotchFreq = FiltCands(1)/Fnyq + 0.0075
 end
 
if minNormNotchFreq <= 0.035 % equiripple highpass filter
[actRp,actAs,b] = LVxDesignEquirippHPF(Rp,As,maxNormNotchFreq,maxNormNotchFreq+0.07);
elseif maxNormNotchFreq >= 0.95 % equiripple lowpass filter
[actRp,actAs,b] = LVDesignEquirippLPF(Rp,As,minNormNotchFreq-0.05,minNormNotchFreq);
else % equiripple notch filter
        wp1 = minNormNotchFreq - 0.035
        ws1 = minNormNotchFreq
        ws2 = maxNormNotchFreq
        wp2 = maxNormNotchFreq + 0.035
[actRp,actAs,b] = LVxDesignEquirippNotch(Rp,As,wp1,ws1,ws2,wp2);
end
 a = 1;   
    % do the filtering in time domain of the entire audio signal
Output = filter(b,a,TstSig);

lenTDConvOutput = length(Output);

Output = Output*(1/max(abs(Output)));
Fr = abs(fft(Output,n));
% do frequency domain filtering using the filter
imp = b;
figure(999)
stem(imp)
xlabel('Sample')
ylabel('Amplitude')

if OvrLap==0
     ProperFDMat = FDMat4Recon;
 else
     ProperFDMat = lgFty;
 end
szlgFty = size(ProperFDMat);

filtfft = fft(imp,2*SzWin);
szfiltfft = size(filtfft);
filtTemp = (conj(filtfft'))*ones(1,szlgFty(2));

szPropFDMat = size(ProperFDMat)
szFiltTemp = size(filtTemp)

fdOutMat = ProperFDMat.*filtTemp;
szfdOutMat = size(fdOutMat);
%==========================================================================
figure(998)
yy = abs(fdOutMat);
yy = yy/max(max(yy));
yy = 20*log10(yy+10^(-6));
szyy = size(yy);
yy = yy(1:fix(szyy(1)/2)+1,:);

mesh(yy)
xlabel('Frame (Filtered Output Signal)')
ylabel('Bin')
zlabel('Mag, dB')
%==========================================================================
lgTD = real(ifft(fdOutMat));
szlgTD = size(lgTD);
newTD = lgTD(1:SzWin,1:szlgTD(2));
newTD(:,2:szlgTD(2)) = newTD(:,2:szlgTD(2)) + lgTD(SzWin+1:szlgTD(1),1:szlgTD(2)-1);
newTD = newTD';
%========================================================================
else  % didn't do any filtering as no candidates were found
    Output = TstSig;
    if OvrLap==0
     ProperFDMat = FDMat4Recon;
 else
     ProperFDMat = lgFty;
    end
szlgFty = size(ProperFDMat);

fdOutMat = ProperFDMat;
szfdOutMat = size(fdOutMat);
lgTD = real(ifft(fdOutMat));
szlgTD = size(lgTD);
newTD = lgTD(1:SzWin,1:szlgTD(2));
newTD = newTD';
end

figure(189)
xvec = (Fs/2)*[0:1:length(Fr)/2]/(length(Fr)/2);
plot(xvec,Fr(1,1:fix(length(Fr)/2)+1))
xlabel('Frequency, Hz (Filtered Entire Signal)')
ylabel('Magnitude')
pause(4)
sound(Output,8000)

% reconstruct signal frames==========================================

tdSig = LVxTDMat2SigVec(newTD,SzWin,SampOverLap);
%==========================================================================
tdSig = tdSig*(1/max(abs(tdSig)));
tdSig  = tdSig(1,1:lenFile);
pause(6)

lenFDConvOutput = length(tdSig);

sound(tdSig,8000)

n = 2^(ceil(log2(length(tdSig))));
finFr = abs(fft(tdSig,n));
figure(190)
xvec = (Fs/2)*[0:1:length(finFr)/2]/(length(finFr)/2);
plot(xvec,finFr(1,1:fix(length(finFr)/2)+1))
xlabel('Frequency, Hz (Filtered Signal, Reconstructed by Frames)')
ylabel('Magnitude')




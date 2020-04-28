function [ToneFreq,Fnyq,BnSp] = LVx_DetnFiltIIRContTone(A,Freq,RorSS,SzWin,OvrLap,AudSig,Rp,FltOrd)
% function [ToneFreq,Fnyq,BnSp] = LVx_DetnFiltIIRContTone(A,Freq,RorSS,SzWin,OvrLap,AudSig,Rp,FltOrd)
%
% This program evaluates the list of candidate interfering tones in the same manner 
% as the script LVx_DetnFiltFIRContTone, with the exception that filtering is accomplished
% using a Chebyshev Type-I filter (lowpass, notch, or highpass, depending on the location of the
% interfering frequency, with maximum passband ripple of Rp dB and order
% FltOrd). The first six input arguments and the output arguments of this
% script are identical to those of LVx_DetnFiltFIRContTone.
% Test calls:
% [T,F,B] = LVx_DetnFiltIIRContTone(0.015,200,0,512,1,1,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.01,210,0,512,1,1,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.01,200,0,512,1,2,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.1,200,0,512,1,3,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.01,500,0,512,1,1,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.01,500,0,512,1,2,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.08,500,0,512,1,3,1,2)
%
% [T,F,B] = LVx_DetnFiltIIRContTone(0.025,200,1,512,1,1,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.02,200,1,512,1,2,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.1,200,1,512,1,3,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.02,500,1,512,1,1,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.018,500,1,512,1,2,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.1,500,1,512,1,3,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.07,150,0,512,1,3,1,2)
%
% [T,F,B] = LVx_DetnFiltIIRContTone(0.005,1000,0,512,0,1,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.008,1000,0,512,1,1,1,2)
% [T,F,B] = LVx_DetnFiltIIRContTone(0.005,1000,0,512,1,1,1,2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
ftsz1 = 16;
ftsz2 = 14;
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
%set(gca,'fontsize',ftsz2)
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
clf
hold on
stem(meanframes)
hline = line([0,length(meanframes)],[xhist(histthresh),xhist(histthresh)]);
%set(hline,'linestyle',':')
xlabel('Frame (Entire Signal)')
ylabel('Avg Bin Mag')

figure(98)
clf
mesh(logfabsfty)
xlabel('Frame (Entire Signal)')
ylabel('Bin')
zlabel('Mag, dB')

figure(99)
clf
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
clf
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
lenFiltCands = length(FiltCands);
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
maxFrqDiff = (FiltCands(adjFrqs(lenadjFrqs)) - FiltCands(adjFrqs(1)));
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
% N = 2;
 if length(FiltCands)>1
    minNormNotchFreq = FiltCands(adjFrqs(1))/Fnyq;
    maxNormNotchFreq = FiltCands(adjFrqs(lenadjFrqs))/Fnyq;
 else
     minNormNotchFreq = FiltCands(1)/Fnyq - 0.025;
     maxNormNotchFreq = FiltCands(1)/Fnyq + 0.025;
 end
    if minNormNotchFreq <= 0.035 % highpass filter
      [b,a] = cheby1(2*FltOrd,Rp,[maxNormNotchFreq + 0.025],'high');
    elseif maxNormNotchFreq >= 0.95 % lowpass filter
         minNormNotchFreq = minNormNotchFreq - 0.025;
     [b,a] = cheby1(2*FltOrd,Rp,minNormNotchFreq);
    else % notch filter
     [b,a] =  cheby1(FltOrd,Rp,[minNormNotchFreq-0.015,maxNormNotchFreq+0.015],'stop');
    end
    % do the filtering in time domain of the entire audio signal
Output = filter(b,a,TstSig);

lenTDConvOutput = length(Output);

Output = Output*(1/max(abs(Output)));
Fr = abs(fft(Output,n));
% do frequency domain filtering using the filter
imp = filter(b,a,[1,zeros(1,2*SzWin-1)]);
figure(999)
stem(imp(1,1:60))
xlabel('Sample (Filter Impulse Response)')
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
% reconstruct original test signal from frames 
fdOutMat = ProperFDMat;
szfdOutMat = size(fdOutMat);
lgTD = real(ifft(fdOutMat));
szlgTD = size(lgTD);
newTD = lgTD(1:SzWin,1:szlgTD(2));
%newTD(:,2:szlgTD(2)) = newTD(:,2:szlgTD(2)) + lgTD(SzWin+1:szlgTD(1),1:szlgTD(2)-1);
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

if ~(DoFilt==0)
figure(777)
FrImp = abs(fft(imp,2^12));
FrImp = 20*log10(FrImp+10^(-8));
xvec = [0:1:(2^11)]/(2^11);
plot(xvec,FrImp(1,1:(2^11 + 1)))
xlabel('Frequency, Units of \pi')
ylabel('Magnitude')
end




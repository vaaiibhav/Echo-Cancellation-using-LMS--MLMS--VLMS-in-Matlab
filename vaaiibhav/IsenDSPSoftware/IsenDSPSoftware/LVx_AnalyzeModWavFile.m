function [ToneFreq] = LVx_AnalyzeModWavFile(strWavFile,SzWin,OvrLap,FrqRg,A,Freq,Rp,As,Act,strOutWavFile)
% [T] =
% LVx_AnalyzeModWavFile('drwatson8Kplus400Hz.wav',512,0,[200,600],0,0,1,20,1,[])
% function [ToneFreq] = LVxAnalyzeModWavFile(strWavFile,SzWin,OvrLap,FrqRg,A,Freq,Rp,As,Act,strOutWavFile)
% strWavFile is a .wav file in the target folder (or a folder on the MathScript search path
% to be read and processed according the value of Act.
% strWavFile is a .wav file in the target folder to be read and processed
% according the value of Act.
% Act = 1 means analyze only and report candidate interfering frequencies as output variable T
% Act = 2 means analyze and filter candidate frequencies if they lie within the frequency range designated
% by FrqRg (see below)
% Act = 3 means analyze, filter, and save the filtered output as a .wav file under the name strOutWavFile
% Act = 4 means do not analyze or filter, but make a .wav file from the test signal, which might have an
% added component if A has magnitude greater than zero. The action may be used to generate .wav files
% having specific frequencies of interfering sinusoid. 
% strOutWavFile must be different from strWavFile
%
% TonFreq, SzWin, OvrLap, A, and Freq have the same meaning and use as in the script
% LVx_DetectContTone, LVx_DetnFiltFIRContTone, and LVx_DetnFiltIIRContTone
%
% FrqRg is a vector of two frequencies in Hz that define a range over which
% to consider filtering. This useful when there are multiple interfering
% frequencies since only simple filters (LPF, HPF, and Notch) can be used
% per call to this script. Successive calls specifying different frequency ranges
% will allow removal, for example, of multiple interfering frequencies
% lying at a distance from one another.
%
% Rp and As are the desired maximum passband ripple in dB and minimum
% stopband attenuation in dB for the filter. Generally, the minimum value of As should be experimentally 
% determined and used; due to masking caused by transient signal tones, the audibility of a low level 
% persistent tone can often be eliminated using As = 20 dB or 30 dB. Very strong persistent tones may 
% require correspondingly higher values of As.
% Test calls:
%
% [T] = LVx_AnalyzeModWavFile('drwatson8Kplus400Hz.wav',512,0,[200,600],0.025,400,1,20,2,[])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
ftsz1 = 16;
ftsz2 = 14;
fr = [];

if ((Act==2|Act==3) & (~(length(FrqRg)==2)))
 error('FrqRg must have two frequencies, in Hz, constituting a frequency range for possible filtering')
end
[y,Fs,bits] = wavread(strWavFile);

NyqRate = Fs/2;
Fnyq = NyqRate;
y = y'; 
lenFile = length(y);
y = (1/max(abs(y)))*y;
t = [0:1:lenFile-1]/Fs;

TstSig = y;
maxabsTstSig = max(abs(TstSig))

if ~(rem(SzWin,2)==0)
    error('SzWin must be an even number! Exiting...')
    return
end
TstSig = TstSig*(0.99/max(abs(TstSig)));
TstSig = TstSig + A*cos(2*pi*t*Freq);
TstSig = TstSig*(0.99/max(abs(TstSig)));
maxabsTstSig = max(abs(TstSig))


if Act==4 % no filtering, make .wav file with TstSig
    if isempty(strOutWavFile)
        error('strOutWavFile is empty')
    end
    TF = strcmpi(strWavFile,strOutWavFile);
    if TF==1
        error('Name of output .wav file is same as name of input .wav file')
    else
    wavwrite(TstSig,Fs,bits,strOutWavFile)
    ToneFreq = [];
    return
    end  
end

n = 2^(ceil(log2(lenFile)))
Fr = abs(fft(TstSig,n));

figure(188)

xvec = (Fs/2)*[0:1:length(Fr)/2]/(length(Fr)/2);
plot(xvec,Fr(1,1:fix(length(Fr)/2)+1))
set(gca,'fontsize',ftsz2)
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

%figure(122)
%clf
%hold on
%stem(meanframes)
%hline = line([0,length(meanframes)],[xhist(histthresh),xhist(histthresh)]);
%set(hline,'linestyle',':')
%xlabel('Frame (Entire Signal)')
%ylabel('Avg Bin Mag')

%figure(98)
%clf
%mesh(logfabsfty)
%xlabel('Frame (Entire Signal)')
%ylabel('Bin')
%zlabel('Mag, dB')

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
%figure(554)
%clf
%xvec = [0:1:length(avgDeriv)-1]/(length(avgDeriv)-1)*(Fs/2);
%stem(xvec,avgDeriv)
%xlabel('Frequency, Hz')
%ylabel('Avg Mag F.O.D. for all Frames')

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

%figure(555)
%clf
%xvec = [0:1:length(normDerivMM)-1]/(length(normDerivMM)-1)*(Fs/2);
%stem(xvec,normDerivMM)
%xlabel('Frequency, Hz')
%ylabel('Avg Mag F.O.D. for low Mag Frames')

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

%NToneFreq = ToneFreq/Fs; % normalized interfering frequencies

sound(TstSig,8000)

if Act==1 % analyze only, don't need code below
return
end
FrqRg = sort(FrqRg);
if Act==2|Act==3 % analyze and filter
% Filter out a selected frequency or group of frequencies.
FiltCands = ToneFreq;
DCFrqs = find(FiltCands <= FrqRg(2) & FiltCands >= FrqRg(1)); % consider only 
% frequencies lying with FrqRg for filtering
FiltCands = FiltCands(DCFrqs);
lenFiltCands = length(FiltCands);
frDiffs = FiltCands(1,2:lenFiltCands) - FiltCands(1,1:lenFiltCands-1);
adjFrqs = find(frDiffs==BnSp);

if ~(isempty(adjFrqs)) 
lenadjFrqs = length(adjFrqs);
adjFrqs = [adjFrqs,(adjFrqs(lenadjFrqs)+1)];
lenadjFrqs = length(adjFrqs)
maxFrqDiff = (FiltCands(adjFrqs(lenadjFrqs)) - FiltCands(adjFrqs(1)))
    if maxFrqDiff <= 50
        DoFilt = 1;
        comment = 'max freq diff <= 50'
    else
        DoFilt = 0;
        comment = 'max freq diff > 50'
    end
else % no adjacent frequencies
    if length(FiltCands)==0
        DoFilt=0;
    else 
        FiltCands = FiltCands(length(FiltCands))
        DoFilt = 1;
    end
end
DoFilt = DoFilt
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
%figure(999)
%stem(imp)
%xlabel('Sample')
%ylabel('Amplitude')

if OvrLap==0
     ProperFDMat = FDMat4Recon;
 else
     ProperFDMat = lgFty;
 end
szlgFty = size(ProperFDMat);

filtfft = fft(imp,2*SzWin);
szfiltfft = size(filtfft);
filtTemp = (conj(filtfft'))*ones(1,szlgFty(2));

szPropFDMat = size(ProperFDMat);
szFiltTemp = size(filtTemp);

fdOutMat = ProperFDMat.*filtTemp;
szfdOutMat = size(fdOutMat)

figure(998)
clf
yy = abs(fdOutMat);
yy = yy/max(max(yy));
yy = 20*log10(yy+10^(-6));
szyy = size(yy);
yy = yy(1:fix(szyy(1)/2)+1,:);

mesh(yy)
xlabel('Frame (Filtered Output Signal)')
ylabel('Bin')
zlabel('Mag, dB')

figure(189)
xvec = (Fs/2)*[0:1:length(Fr)/2]/(length(Fr)/2);
plot(xvec,Fr(1,1:fix(length(Fr)/2)+1))
xlabel('Frequency, Hz (Filtered Entire Signal)')
ylabel('Magnitude')
pause(4)
sound(Output,8000)

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
szfdOutMat = size(fdOutMat)

figure(998)
clf
yy = abs(fdOutMat);
yy = yy/max(max(yy));
yy = 20*log10(yy+10^(-6));
szyy = size(yy);
yy = yy(1:fix(szyy(1)/2)+1,:);
mesh(yy)
xlabel('Frame (Unfiltered Output Signal)')
ylabel('Bin')
zlabel('Mag, dB')
end % for DoFilt



end % (for Act=2 or Act=3 analyze & filter)

if Act==3
    if isempty(strOutWavFile)
        error('strOutWavFile is empty')
    end
    TF = strcmpi(strWavFile,strOutWavFile);
    if TF==1
        error('Name of output .wav file is same as name of input .wav file')
    else
      Output = 0.99*Output/max(abs(Output));
      wavwrite(Output,Fs,bits,strOutWavFile)
    end
end

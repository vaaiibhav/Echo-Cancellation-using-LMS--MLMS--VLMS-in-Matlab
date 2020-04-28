function LVxReverbDereverb(Delay,DecayRate,Mu,PeakSep)
% funtion LVxReverbDereverb(Delay,DecayRate,Mu,PeakSep)
% Delay is the number of samples of delay used in the single feedback stage
% to generate the reverberated sound file.
% DecayRate is a real number which is feedback gain.
% Mu is the usual LMS update weight term.
% PeakSep the the number of samples by which to separate detected peaks in
% the autocorrelation sequence of the reverberated sound.
% Uses a Dual-H architecture to estimate DecayRate, and uses the best estimate of DecayRate
% from the Dual-H process along with the estimated value of Delay (from autocorrelation) 
% in an FIR to deconvolve the entire reverberative signal.
% Does not archive or plot the values of RRLE or current best estimate of
% DecayRate. Plays the reverberative sound once and the dereverberated
% sound once. These sounds may be played again after running the script by
% making the calls global ReverbSnd and global DereverbSnd in the Command
% window and then making the calls sound(ReverbSnd,8000) and
% sound(DereverbSnd,8000) at will. Uses the audio file 'drwatsonSR4K.wav' 
% as the test signal.   
% Sample calls: 
% LVxReverbDereverb(325,0.7,0.05,200)
% LVxReverbDereverb(1800,0.975,0.05,200)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

global ReverbSnd
global DereverbSnd

ReverbSnd = [];
DereverbSnd = [];

Mu = real(Mu);
PeakSep = fix(abs(real(PeakSep)));
Delay = abs(real(Delay));
if Delay>4600
   Comment = 'Delay is limited to 4600; setting Delay to 4600'
Delay = 4600;
end

DecayRate = real(DecayRate);
if abs(DecayRate)>1.0 
   Comment = 'Magnitude of DecayRate must be not greater than 1.0'
   return
end

[y,Fs,bits] = wavread('drwatsonSR4K.wav');

y = y';
lenOrigFile = length(y);

lenFile = length(y);
y = (1/max(abs(y)))*y;

NetDrWatsonDurInSecs = 2^17/Fs;
BinWidth = 1/NetDrWatsonDurInSecs;

ReverbSnd = y(1,1:lenFile);
ReverbSnd = filter(1,[1 zeros(1,Delay-1) -DecayRate],ReverbSnd);
ReverbSnd = (0.99/max(abs(ReverbSnd)))*ReverbSnd;
DereverbSnd = zeros(1,length(ReverbSnd));
%==================================================================================

aaa = abs(fft(y,2^17));
bbb = abs(fft(ReverbSnd,2^17));

figure(3957)

limfreqplot = 4500;

subplot(311)
plotlim1 = 1.2*max(abs(aaa(1,1:limfreqplot)));
plot( BinWidth*(1:1:limfreqplot), aaa(1,1:limfreqplot))
ylabel(['Magnitude'])
xlabel(['(a) Frequency (Hz), Orig Audio Signal'])
axis([0, BinWidth*limfreqplot, 0, plotlim1])

subplot(312)
plotlim2 = 1.2*max(abs(bbb(1,1:limfreqplot)));
plot(BinWidth*(1:1:limfreqplot),bbb(1,1:limfreqplot))
ylabel(['Magnitude'])
if DecayRate>0
xlabel(['(b) Freq (Hz), Orig Signal Conv w/ Xfer Fcn 1/(1 - ',num2str(DecayRate),'e^{-j\omega',num2str(Delay),'})'])
else
xlabel(['(b) Freq (Hz), Orig Signal Conv w/ Xfer Fcn 1/(1 + ',num2str(abs(DecayRate)),'e^{-j\omega',num2str(Delay),'})'])
end
axis([0, BinWidth*limfreqplot, 0, plotlim2])

freqvec = 0:pi/((2^16)-1):pi;
lenFreqVec = length(freqvec);
freqResp = abs(ones(1,lenFreqVec)./(ones(1,lenFreqVec) - DecayRate*exp(-j*Delay*freqvec)));

subplot(313)
plotlim3 = 1.2*max(abs(freqResp(1,1:limfreqplot)));
plot(BinWidth*(1:1:limfreqplot),freqResp(1,1:limfreqplot))
ylabel(['Magnitude'])
if DecayRate>0
xlabel(['(c) Freq Resp (Hz), Xfer Fcn 1/(1 - ',num2str(DecayRate),'e^{-j\omega',num2str(Delay),'})'])
else
xlabel(['(c) Freq Resp (Hz), Xfer Fcn 1/(1 + ',num2str(abs(DecayRate)),'e^{-j\omega',num2str(Delay),'})'])
end
axis([0, BinWidth*limfreqplot, 0, plotlim3])

figure(771);

subplot(311)
plot(ReverbSnd,'b')
ylabel(['Amplitude'])
xlabel(['(a) Sample, Conv Signal; Delay = ',...
      num2str(Delay/(Fs)),' Sec.'])
axis([0  length(ReverbSnd) -1 1])

b = xcorr(ReverbSnd,fix(1.5*Delay)); % Compute the autocorrelation sequence
bAbs = abs(b);
b = (0.99/max(abs(b)))*b;
[absloc, absval] = LVfindPeaks(bAbs,2,PeakSep); 

estDelay = abs(absloc(1) - absloc(2)) 
%===================LMS Filtering to get DecayRate================================
mult = 1;
lenFile = fix(mult*length(ReverbSnd));
BestTapWt(1:lenFile) = 0;
BestTapWt(1:estDelay +1) = 0;
TestErr(1:lenFile) = 0;
CurBestRRLE(1:lenFile) = ones(1,lenFile)*10^-6;
CurBestRRLE(1:estDelay+11) = 1;
TestRRLE = 1;
NewTestTapWt = 0;
Ctr = estDelay+11;

while Ctr<lenFile & log10(CurBestRRLE(Ctr-1))<5
if TestRRLE > CurBestRRLE(Ctr-1)
   BestTapWt(Ctr) = NewTestTapWt;
   CurBestRRLE(Ctr) = TestRRLE;
else
   CurBestRRLE(Ctr) = CurBestRRLE(Ctr-1);
   BestTapWt(Ctr) = BestTapWt(Ctr-1);
end     
TestFiltOut = NewTestTapWt*ReverbSnd(Ctr-estDelay);      	  
TestErr(Ctr) = ReverbSnd(Ctr) - TestFiltOut;
TestRRLE = sum(ReverbSnd(Ctr-10:Ctr).^2)/(sum(TestErr(Ctr-10:Ctr).^2)+10^-6);
NewTestTapWt = NewTestTapWt + Mu*TestErr(Ctr)*ReverbSnd(Ctr-estDelay)/(sum(ReverbSnd(Ctr-estDelay-10:Ctr-estDelay).^2)+10^-6);  
Ctr = Ctr + 1;
end
[theLocs, thePeaks] = findPeaks(CurBestRRLE,1,1);
estDecayRate = BestTapWt(theLocs(1))

figure(80)

subplot(2,1,1)
plot(BestTapWt(1,1:Ctr-1),'b');
hold on
hndPlot = plot(Ctr-1,BestTapWt(Ctr-1));
grid on
hold off
xlabel(['(a) Sample, Tap Wt (DecayRate); Est = ',num2str(estDecayRate)])
ylabel(['Amplitude'])
axis([0 Ctr-1 (min(BestTapWt) - 0.2) (max(BestTapWt)+0.2)])

subplot(2,1,2)
theans = 20*log10(CurBestRRLE(1,1:Ctr-1)+10^-5);
plot(theans)
hold on
plot(Ctr-1,theans(Ctr-1));
grid on
ylabel(['Mag, dB'])
xlabel(['(b) Sample, RRLE. (Max(RRLE) = ',num2str(20*log10(thePeaks(1))),' dB)'])
axis([0,Ctr-1,-5,120])
hold off
% ============================================================

DereverbSnd = filter([1, zeros(1,estDelay-1),  -estDecayRate],1,ReverbSnd); % easier way
DereverbSnd = (0.99/max(abs(DereverbSnd)))*DereverbSnd;

figure(771)
subplot(313)
plot(DereverbSnd,'b')
ylabel(['Amplitude'])
xlabel(['(c) Sample, Deconvolved Audio Signal'])
axis([0  length(DereverbSnd)  -1 1])

subplot(312)
plot(b,'b')
ylabel(['Amplitude'])
xlabel(['(b) AutoCorr Seq; Est/Delay = ',num2str(estDelay/Fs),' Sec.'])
axis([0  length(b)  -1.1  1.1])

figure(771)

ReverbSnd = interp(ReverbSnd,2,10,0.5);
sound(ReverbSnd,8000);
pause(5)
DereverbSnd = interp(DereverbSnd,2,10,0.5);
sound(DereverbSnd,8000);

Comment = 'global sound output variable names are ReverbSnd and DereverbSnd'
Comment = ['Output sound sample rate = ',num2str(2*Fs)]

return

function [Locs, Pks] = LVfindPeaks(Mat,QuanPks,PeakSep)
% function [Locs, Pks] = LVfindPeaks(Mat,QuanPks,PeakSep)
% Mat is a vector or matrix for which QuanPks peak values are sought, each 
% of which is at least PeakSep samples distant from the next closest peak 
% value. Locs tells the locations of the peak values in the matrix. To
% make adjacent samples of Mat eligible to be peaks, use PeakSep = 0.
% Sample call:
% [Locs, Pks] = LVfindPeaks([0:1:7],8,0)
% [Locs, Pks] = LVfindPeaks([0:1:7],8,1)
% [Locs, Pks] = LVfindPeaks((cos(2*pi*3*[0:1:63]/64)),4,8)

aa = size(Mat); 
Mat = Mat(:);  
nums = find(isfinite(Mat));
minVal = min(Mat)-1;
notnums = find(isnan(Mat));
Mat(notnums) = minVal;
    if QuanPks > length(nums), 
        QuanPks = length(nums);
     Comment = 'Num peaks reduced'; 
    end
Pks = zeros(1,QuanPks);
Locs = zeros(1,QuanPks);

for ctr = 1:1:QuanPks 
  [aa,bb] =  max(Mat);
  Pks(1,ctr) = aa;
  Locs(1,ctr) = bb;
  minLocToMin = max([1, (Locs(ctr)-PeakSep)]);
  maxLocToMin = min([length(Mat),(Locs(ctr)+ PeakSep)]);
  locsToMin = minLocToMin:1:maxLocToMin;
  Mat(locsToMin) = minVal; 
end;

NumNotPeaks = length(find(Pks==minVal));
Pks = Pks(1,1:length(Pks)-NumNotPeaks);
Locs = Locs(1,1:length(Locs)-NumNotPeaks);


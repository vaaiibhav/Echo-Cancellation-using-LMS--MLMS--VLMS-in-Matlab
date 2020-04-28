function LVxLMSvNLMS(kVec,NoTaps,Delta,Mu,SigType,Freq,Delay)
% function LVxLMSvNLMS(kVec,NoTaps,Delta,Mu,SigType,Freq,Delay)
%
% kVec is a vector of amplitudes to impose on the test signal as a succession
% of equally spaced amplitude segments over a test signal of 15 times the 
% filter length.
% NoTaps is the number of taps to use in the LMS adaptive filter;
% The value for Stepsize is computed as Delta/(P + sm) where
% P = sum(testSig.^2) for the samples of testSig in the filter, and this
% is used in the NLMS algorithm (sm is a small number such as 10^(-6)
% Simultaneously, the test signal is processed by an LMS filter of the same
% length using a constant stepsize equal to Mu.
% SigType 0 gives a DC test signal,
% SigType 1 = random noise, 2 = cosine of frequency Freq
% SigType 3 = Nyquist rate, and 4 = Half-Band Frequency
% The current filter output sample is subtracted from the current
% test signal sample, delayed by Delay samples. Delay must be between
% 1 and NoTaps for the error signal to be able to go to zero
% Three plots are created, one of the test signal, one of the 
% LMS adaptive FIR error/output signal, and one of the NLMS error signal.
% A typical call:
% LVxLMSvNLMS([1,5,15,1],100,0.5,0.00005,1,[15],50)
% LVxLMSvNLMS([1,5,25,1],100,0.5,0.00005,1,[15],50)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

kVec = real(kVec);
NoTaps = floor(abs(NoTaps));
Err = [];

if NoTaps > 500
   NoTaps = 500;
end

SR = 15*NoTaps;
t = 0:1/SR:1-1/SR;

AmpProfile = ones(1,SR);
DivPts = [1,fix(SR*[1:1:(length(kVec)-1)]/length(kVec))];
DivPts = [DivPts,(SR+1)];

for Ctr = 1:1:length(kVec)
AmpProfile(1,DivPts(Ctr):DivPts(Ctr+1)-1) = kVec(Ctr);
end

%AmpProfile;

if SigType==0
DataVec = AmpProfile.*ones(1,SR);
elseif SigType==1
DataVec = AmpProfile.*randn(1,SR);
elseif SigType==2
DataVec = AmpProfile.*cos(2*pi*t*Freq);
elseif SigType==3
DataVec = AmpProfile.*(-1).^(0:1:SR-1);
elseif SigType==4
DataVec = AmpProfile.*cos(2*pi*t*(length(t)/4));    
else
end

TapWt(1:NoTaps,1) = 0;
NTapWt(1:NoTaps,1) = 0;
FiltOut(1,1:length(DataVec)) = 0;
NFiltOut(1,1:length(DataVec)) = 0;
Err(1,1:length(DataVec)) = 0;
NErr(1,1:length(DataVec)) = 0;
n = 1:1:NoTaps;
 
for CurDtaPtr = NoTaps+1:1:length(DataVec)
% NLMS
NFiltOut(1,CurDtaPtr) = DataVec(CurDtaPtr-n)*NTapWt(n,1);  
NErr(1,CurDtaPtr)= DataVec(1,CurDtaPtr-Delay) - NFiltOut(1,CurDtaPtr); 
PwrNFilt = sum(DataVec(CurDtaPtr-n).^2) + 10^(-6);
NetDelta = Delta/PwrNFilt;
NTapWt(n,1)= NTapWt(n,1)+ NetDelta*NErr(1,CurDtaPtr)*DataVec(1,CurDtaPtr-n)';
% LMS below
FiltOut(1,CurDtaPtr) = DataVec(CurDtaPtr-n)*TapWt(n,1);  
Err(1,CurDtaPtr)= DataVec(1,CurDtaPtr-Delay) - FiltOut(1,CurDtaPtr); 
TapWt(n,1)= TapWt(n,1)+ Mu*Err(1,CurDtaPtr)*DataVec(1,CurDtaPtr-n)';
end

figure(17)

subplot(321)
rDataVec = DataVec;
plot(rDataVec)
zzc = max(abs(DataVec));
ylabel(['Amplitude'])
xlabel(['(a) Sample, Test Signal'])
axis([0, inf, min([0,1.2*min(rDataVec)]), 1.2*max(rDataVec)])

subplot(322)
rDataVec = DataVec;
plot(rDataVec)
zzc = max(abs(DataVec));
ylabel(['Amplitude'])
xlabel(['(b) Sample, Test Signal'])
axis([0, inf, min([0,1.2*min(rDataVec)]), 1.2*max(rDataVec)])


subplot(323)

plot(Err)
vcx = max(abs(Err));
ylabel(['Amplitude'])
xlabel(['(c) Sample (Output/Error, LMS)'])
axis([0, SR, -1.2*max(Err), 1.2*max(Err)])

subplot(324)
plot(NErr)
vcx = max(abs(NErr));
ylabel(['Amplitude'])
xlabel(['(d) Sample (Output/Error, NLMS)'])
axis([0, SR, -1.2*max(NErr), 1.2*max(NErr)])

subplot(325)
stem(TapWt,'bo')
ylabel(['Amplitude'])
xlabel(['(e) Tap No. (LMS)'])
axis([0, NoTaps, -1.2*max(TapWt), 1.2*max(TapWt)])

subplot(326)
stem(NTapWt,'bo')
ylabel(['Amplitude'])
xlabel(['(f) Tap No. (NLMS)'])
axis([0, NoTaps, -1.2*max(NTapWt), 1.2*max(NTapWt)])







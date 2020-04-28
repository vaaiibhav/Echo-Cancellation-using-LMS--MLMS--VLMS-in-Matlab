function LVxLMSAdaptFiltMuCheck(k,NoTaps,C,SigType,Freq,Delay)
% function LVxLMSAdaptFiltMuCheck(k,NoTaps,C,SigType,Freq,Delay)
%
% k is an amount of white noise to use;
% NoTaps is the number of taps to use in the LMS adaptive filter;
% The value for Stepsize is computed as C/(10*P*NoTaps) where
% P = sum(testSig.^2)/length(testSig);
% SigType 0 gives a DC test signal,
% SigType 1 = random noise, 2 = cosine of frequency Freq
% SigType 3 = Nyquist rate, and 4 = Half-Band Frequency
% The current filter output sample is subtracted from the current
% test signal sample, delayed by Delay samples. Delay must be between
% 1 and NoTaps for the error signal to be able to go to zero
% Two plots are created, one of the test signal, and one of the 
% LMS adaptive FIR error/output signal
% A typical call:
% LVxLMSAdaptFiltMuCheck(1,100,1,1,[15],50)
% LVxLMSAdaptFiltMuCheck(1,100,4,1,[15],50)
% LVxLMSAdaptFiltMuCheck(1,100,6,1,[15],50)
% LVxLMSAdaptFiltMuCheck(1,100,8,1,[15],50)
% LVxLMSAdaptFiltMuCheck(1,100,9.5,1,[15],50)
% LVxLMSAdaptFiltMuCheck(1,100,10,1,[15],50)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

k = real(k);
NoTaps = floor(abs(NoTaps));
Err = [];

if NoTaps > 500
   NoTaps = 500;
end

SR = 15*NoTaps;
t = 0:1/SR:1-1/SR;

if SigType==0
DataVec = k*ones(1,SR);
elseif SigType==1
DataVec = k*randn(1,SR);
elseif SigType==2
DataVec = k*cos(2*pi*t*Freq);
elseif SigType==3
DataVec = k*(-1).^(0:1:SR-1);
elseif SigType==4
DataVec = k*cos(2*pi*t*(length(t)/4));    
else
end

P = (1/length(DataVec))*sum(DataVec.^2);
Mu = C/(10*NoTaps*P);

TapWt(1:NoTaps,1) = 0;
FiltOut(1,1:length(DataVec)) = 0;
Err(1,1:length(DataVec)) = 0;
n = 1:1:NoTaps;
 
for CurDtaPtr = NoTaps+1:1:length(DataVec)
FiltOut(1,CurDtaPtr) = DataVec(CurDtaPtr-n)*TapWt(n,1);  
Err(1,CurDtaPtr)= DataVec(1,CurDtaPtr-Delay) - FiltOut(1,CurDtaPtr); 
TapWt(n,1)= TapWt(n,1)+ 2*Mu*Err(1,CurDtaPtr)*DataVec(1,CurDtaPtr-n)';
end

figure(17)

subplot(311)
stem(DataVec(1,1:200))
zzc = max(abs(DataVec));
ylabel(['Amplitude'])
xlabel(['(a) Sample (First 200 Samples of Test Signal)'])
axis([0, inf, -1.2*max(DataVec), 1.2*max(DataVec)])

subplot(312)

plot(Err)
vcx = max(abs(Err));
ylabel(['Amplitude'])
xlabel(['(b) Sample (Output/Error)'])
axis([0, SR, -1.2*max(Err), 1.2*max(Err)])

subplot(313)
stem(TapWt,'bo')
ylabel(['Amplitude'])
xlabel(['(c) Tap Number'])
axis([0, NoTaps, -1.2*max(TapWt), 1.2*max(TapWt)])





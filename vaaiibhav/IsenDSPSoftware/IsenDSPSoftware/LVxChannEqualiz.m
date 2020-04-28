function LVxChannEqualiz(ChImp,Mu,tSig,lenLMS)
% function LVxChannEqualiz(ChImp,Mu,tSig,lenLMS)
% ChImp is the impulse response of the channel which is to be equalized by
% the LMS FIR;
% Mu is the usual update weight
% tSig designates the type of test signal to be used: 0 = random noise
% 1 = Unit impulse sequence of length lenLMS, repeated at least ten times,
% with alternating sign, +, -, +, etc
% 2 = Unit impulse sequence of length lenLMS, repeated at least ten times
% 3 = Sinusoid at the half-band frequency
% lenLMS is the length of the adaptive LMS FIR
% Example calls:
% LVxChannEqualiz([-0.05,0,-0.5,0,-1,0,1,0,0.5,0,0.05],1.5,2,131) % bandpass
% LVxChannEqualiz([1,0,1],2.1,2,91) % notch
% LVxChannEqualiz([1,0,1],2.4,2,191) % notch
% LVxChannEqualiz([1 1],2.2,2,65) % lowpass
% LVxChannEqualiz([1,0,-1],1.8,2,85) % bandpass
% LVxChannEqualiz([1,-1],2.5,2,179) % highpass
% LVxChannEqualiz([1,-0.8],2.3,2,35) % highpass
% LVxChannEqualiz([1,0.8],2.4,2,35) % lowpass
% LVxChannEqualiz([0.3,1,0.3],2.3,2,9) % lowpass
% LVxChannEqualiz([1,0.7],2.3,2,25) % lowpass
% LVxChannEqualiz([1,0,0,0,0,0,1],1.6,2,179) % comb filter
% LVxChannEqualiz([1,0,0,-0.3,0,0,0.05],2.5,2,23) 
% LVxChannEqualiz([1,0,-1],1.5,2,55)
% LVxChannEqualiz( ([1,0,1,0,1,0,1].*hamming(7)'),2.1,2,91)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

figure(1323);

SR = 4096; % for FFT's
lentSig = 15*lenLMS; 
typeWinChann = 'hamming';

ChFiltImpResp = ChImp;
%ChFiltImpResp = (1/max(abs(ChFiltImpResp)))*ChFiltImpResp;
LenChImpResp = length(ChFiltImpResp); 

if tSig==0 % 0 for all noise, 1 for sine waves plus noise
testSig = randn(1,lentSig);
elseif tSig==1
    
testSig = zeros(1,lenLMS);
testSig(1) = 1;
testSig2 = zeros(1,lenLMS);
testSig2(1) = -1;
testSig = [testSig, testSig2];
testSig = [testSig testSig testSig testSig ];   
elseif tSig==2
    testSig = zeros(1,lenLMS);
    testSig(1) = 1;
    testSig = (testSig')*ones(1,4);
    testSig = testSig(:)';    
elseif tSig==3
    t = 0:1/lentSig:1-1/lentSig;
    testSig = sin(2*pi*t*fix(lentSig/4)); 
else
    testSig = randn(1,lentSig);
end

while (length(testSig) < lentSig)
testSig = [testSig testSig];   
end
L = length(testSig);

L = length(testSig);

% ======================================================================================
rawChannelSignal = testSig;
ChannelSignal = (1/max(abs(rawChannelSignal)))*rawChannelSignal;
LenChSig = length(ChannelSignal);
%=======================================================================================

subplot(231)
plot(ChannelSignal,'k')
thisplotlim = 1.2*max(abs(ChannelSignal));
xlabel(['(a) Test Signal'])
ylabel('Amp')
axis([0 length(ChannelSignal) -thisplotlim thisplotlim])
%=======================================================================================
ChannelSignal = conv(ChFiltImpResp,ChannelSignal); % full spectrum shaped by filtering
ChannelSignal = ChannelSignal(1,1:LenChSig);% cut down to same size as desired channel response
ChannelSignal = (1/max(abs(ChannelSignal)))*ChannelSignal;
%==end test signal generation code======================================================

FilterL = lenLMS; 
LenChFiltImpResp = length(ChFiltImpResp);

TapWt = zeros(FilterL,lentSig);
FiltOut = zeros(1,lentSig );
Err = zeros(1,lentSig );

n = 1:1:FilterL; % n is the tap index number.

LenSpecWin = fix(length(ChannelSignal)/10);
LoopStart = FilterL + LenChFiltImpResp+1
LoopLim = lentSig % - LenIdealResp;

%=============================================================================
subplot(232)
stem(ChFiltImpResp,'k');
xlabel('(b) Ch Imp Resp')
ylabel('Amp')
axis([0, LenChFiltImpResp+1, min([0,min(ChFiltImpResp)]),1.1*max(ChFiltImpResp)])

subplot(233)

y = abs(fft(ChFiltImpResp,SR));
y = y(1,1:SR/2+1);
%y = y/max(y);
xvec = 2*[0:1:SR/2]/SR;
plot(xvec,y,'k');
xlabel('(c) Freq, Ch Imp')
ylabel('Mag')
axis([0,1,0,1.1*max(y)])

for CurDtaPtr = LoopStart:1:LoopLim  % This it, the actual computation...   
FiltOut(1,CurDtaPtr) = sum((TapWt(n,CurDtaPtr)').*ChannelSignal(CurDtaPtr-n+1));
% Err(1,CurDtaPtr)=rawChannelSignal(1,CurDtaPtr-fix((FilterL+LenChFiltImpResp)/2))-FiltOut(1,CurDtaPtr);
Err(1,CurDtaPtr)= rawChannelSignal(1,CurDtaPtr-fix(LenChFiltImpResp/2))-FiltOut(1,CurDtaPtr);
%Err(1,CurDtaPtr)= rawChannelSignal(1,CurDtaPtr-(FilterL+LenChFiltImpResp-2)) - FiltOut(1,CurDtaPtr);%=======
%Err(1,CurDtaPtr)= rawChannelSignal(1,CurDtaPtr-(FilterL + LenChFiltImpResp)) - FiltOut(1,CurDtaPtr);%=======
%Err(1,CurDtaPtr)= rawChannelSignal(1,CurDtaPtr-LoopStart+3) - FiltOut(1,CurDtaPtr);%=======

LIMAdjust = 3*(sum(ChannelSignal(1,CurDtaPtr-(FilterL)+1:CurDtaPtr).^2)+0.01);
m = 1:1:FilterL;
TapWt(m,CurDtaPtr+1)= TapWt(m,CurDtaPtr)+ ...
   (2*Mu*Err(1,CurDtaPtr)*((ChannelSignal(1,CurDtaPtr-m+1)))')/LIMAdjust;  
end % end tap update and partition by equalizer/filter design mode============
%============================================================================= 
TapMax = max(max(abs(TapWt)));
  
%=============================================================================   
plotlim3 = 1.2*max(max(abs(TapWt(1:FilterL,CurDtaPtr-1))));
   
subplot(234)
stem(TapWt(1:FilterL,CurDtaPtr-1),'k');
xlabel(['(d) Tap Wts'])
ylabel('Amp')
axis(  [0,FilterL+1,min([-0.25,-plotlim3]),max([0.25,plotlim3])])
%==========================================================================
  
subplot(236)
netImp = conv([ChFiltImpResp],[TapWt(1:FilterL,LoopLim-1)']);
out2 = abs(fft( netImp+10^-12, SR ) );
out2 = out2(1,1:fix(length(out2)/2));
out2 = (1/max(out2))*out2;
lenOut2O2 = length(out2);
xvec = [0:1:lenOut2O2-1]/lenOut2O2;
plot(xvec,out2,'k');
xlabel(['(f) Freq, Equal. Chan'])
ylabel('Mag')
axis([0,1,0,1.1])

SR = 4096;

subplot(235)
out2 = abs(fft(TapWt(1:FilterL,LoopLim-1)',SR ) );
out2 = out2(1,1:SR/2);
lenOut2O2 = length(out2);
xvec = [0:1:lenOut2O2-1]/lenOut2O2;
plot(xvec,out2,'k');
xlabel(['(e) Freq, Tap Wts'])
ylabel('Mag')
axis([0,1,0,1.1*max(out2)])










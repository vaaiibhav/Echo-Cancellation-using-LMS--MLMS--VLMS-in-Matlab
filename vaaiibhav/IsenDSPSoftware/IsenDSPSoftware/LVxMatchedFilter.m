function LVxMatchedFilter(NoiseAmp,TstSeqLen,FlipImpResp)
% function LVxTestReconSineVariablePhase(k1,N,PhaseDegrees)
% Forms a test impulse response, a chirp having a length equal to half of
% TstSeqLen, then builds a test sequence having length TstSeqLen and
% containing the chirp and noise of amplitude NoiseAmp. The test signal is
% then convolved with either the chirp or a time-reversed version of the
% chirp, and the results plotted to demonstrate the principle of matched
% filtering. 
% FlipImpResp: Use 1 for Time-reversed, 0 for not time-reversed impulse
% response.
% Test calls:
% LVxMatchedFilter(0.5,128,0) % imp resp not time reversed
% LVxMatchedFilter(0.5,128,1) % imp resp time reversed
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

SR = 16;
TheLim = floor(1.5*SR);
ImpRespMaxAmp = 0.5;

figure(19997);
clf

% axis([-SR SR -1.2*ImpRespMaxAmp  1.2*ImpRespMaxAmp]) % upper

SR = TstSeqLen;
t = 0:1/(SR/2):1-1/(SR/2);
TheLim = floor(1.5*SR);

stTheType = 'Time-Reversed Chirp';
%-----------------------------------------------------------------
yshift(1,1:1:floor(2.5*TheLim)) = zeros(1,floor(2.5*TheLim));
ConvolSeqVal(1,1:floor(1.5*TheLim)) = zeros(1,floor(1.5*TheLim));
x(1,1:floor(2.5*TheLim)) = zeros(1,floor(2.5*TheLim));
y(1,1:floor(2.5*TheLim)) = zeros(1,floor(2.5*TheLim));

ChirpSR = SR/2;
b = chirp(t,0,1-1/ChirpSR,ChirpSR/2);

timerev = fliplr(b);

if FlipImpResp==1
b = timerev;
end

b = b/(max(abs(b)));
LenImpResp = length(b);
ImpRespMaxAmp = max(abs(b));
x(1,TheLim+SR+1:(TheLim+SR+length(b)) ) = b; %==============================================
xabc = 1:1:floor(2.5*TheLim);
xxabc = 1:1:floor(1.5*TheLim);
y(1,TheLim+1:TheLim+SR/2) = fliplr(chirp(t,0,1-1/ChirpSR,ChirpSR/2));
y = y + NoiseAmp*randn(1,length(y));
ImpRespMaxAmp = max([ abs(b)   abs(y)  ] );

theconv = conv(b,fliplr(y));

subplot(211)

hold off
if SR<65
   xthing = xabc-(TheLim+SR+1);
   stem(xthing,x,'rd'); 
else
plot(xabc-(TheLim+SR+1),x,'r')    
end
hold on

a = 1;
lenY = length(y);
Ysxvec = -(lenY-1):1:0;
if SR<65
stem(Ysxvec,y,'bo'); 
else
plot(Ysxvec,y,'b')    
end

xlabel(['(a) ImpResp (Samp 0 to ',num2str((SR/2) - 1),'), TestSig (Samp Time ',num2str(a-1),')'])
ylabel(['Amplitude'])
axis([-3*SR SR -1.2*ImpRespMaxAmp 1.2*ImpRespMaxAmp])
%==========================================================================

subplot(212)
plot(theconv)
xlabel('(b) Sample')
ylabel('Amplitude')





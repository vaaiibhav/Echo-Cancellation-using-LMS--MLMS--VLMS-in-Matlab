function LVxFreqRespViaZxform(Imp,SR)
% function DemoFreqRespViaZxform(Imp,SR)
%
% Imp is an impulse response in row vector form;
% SR is the desired number of samples for both DTFT and z-transform
% evaluations of the frequency content of the impulse response.
%
% Sample calls: 
%
% LVxFreqRespViaZxform([1 1 1 1 1 1 1 1],512)
% LVxFreqRespViaZxform([1 1 1 1 1 1 1 1],512)
% LVxFreqRespViaZxform([1 0.707 1],512)
% LVxFreqRespViaZxform([0.1 0.2 0.3 0.4 0.4 0.3 0.2 0.1],512)
% LVxFreqRespViaZxform([0.4 0.3 0.2 0.1 0.1 0.2 0.3 0.4],512
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
mag = zeros(1,SR);
TheAng = zeros(1,SR);

for A  = 1:1:SR+1
   z = exp(-j*(A-1)*pi/SR);
   zvec = [z.^(0:1:length(Imp)-1)];  
   VecProd =  sum(Imp.*zvec);  
   mag(A) = abs(VecProd);
   TheAng(A) = angle(VecProd);
end
interAns = fft(Imp,2*SR);
otherAns = abs(interAns);
TheAng = (TheAng);

figure(5679)
clf

subplot(221)
plot((1/SR)*([0:1:SR]),mag)
ylabel(['Magnitude'])
xlabel(['(a) Freq Resp Via z-Xform; X-Axis = NormFreq'])
axis([0 1 0 1.2*max(abs(mag))])
   
subplot(223)
plot((1/SR)*([0:1:SR]),otherAns(1,1:SR+1))  
ylabel(['Magnitude'])
xlabel(['(c) FreqResp Via DTFT; X-Axis = NormFreq'])
axis([0  1  0 1.2*max(abs(otherAns))])

plotlim3up = abs(max(TheAng));
plotlim3down = abs(min(TheAng));

subplot(222)
plot((1/SR)*([0:1:SR]),TheAng)
ylabel(['Radians'])
xlabel(['(b) Phase Resp Via z-Xform; X-Axis = NormFreq'])
axis([0 1  -1.2*plotlim3down 1.2*plotlim3up])

ftphase = (angle(interAns(1,1:SR+1)));
plotlim4up = abs(max(ftphase));
plotlim4down = abs(min(ftphase));

subplot(224)
plot((1/SR)*([0:1:SR]),ftphase)
ylabel(['Radians'])
xlabel(['(d) Phase Resp Via DTFT; X-Axis = NormFreq'])
axis([0 1  -1.2*plotlim4down 1.2*plotlim4up])


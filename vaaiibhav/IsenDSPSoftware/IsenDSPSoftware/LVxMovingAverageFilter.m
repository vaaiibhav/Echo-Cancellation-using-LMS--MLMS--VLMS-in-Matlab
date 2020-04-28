function LVxMovingAverageFilter(MALength,NoiseAmp)
% function LVxMovingAverageFilter(MALength,NoiseAmp)
% MALength is the length of the Moving Average filter
% NoiseAmp is the StdDev of white noise mixed with the coherent
% signal, which consists of five rectangular pulses over 1024 samples,
% each coherent pulse having an amplitude of 1.0 and a
% width of 20 samples.
%
% A typical call:
%
% LVxMovingAverageFilter(20,0.8)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if MALength<1
   MALength = 2;
end

TDSig = [zeros(1,40),ones(1,20),zeros(1,180)];
TDMat = (TDSig')*[ones(1,5)];
TDSig = TDMat(:)';

TDSig = TDSig + NoiseAmp*randn(1,length(TDSig));

x = (1/MALength)*ones(1,MALength);

figure(500)


subplot(221)
stem(x,'bo');
ylabel(['Amplitude'])
xlabel(['(a) Sample'])
axis([0, fix(length(x)+1), -inf, 1.2*max(abs(x))])
   
subplot(223)

plot(TDSig,'b')
ylabel(['Amplitude'])
xlabel(['(c)  Sample'])
axis([0, inf, -1.2*abs(min(TDSig)), 1.2*abs(max(TDSig))])

subplot(224)

TDConv = (1/length(x))*conv(x,TDSig);
plot(TDConv,'b')
ylabel(['Amplitude'])
xlabel(['(d)  Sample'])
axis([0, length(TDConv), -1.2*abs(min(TDConv)), 1.2*abs(max(TDConv))])

subplot(222)

fr = abs(fft(x,4096));
plot([0:1:2048]/2048, fr(1,1:2049))
ylabel(['Magnitude'])
xlabel(['(b)  Freq, Units of \pi'])
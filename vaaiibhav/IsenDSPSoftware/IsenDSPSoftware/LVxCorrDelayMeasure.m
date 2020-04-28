function LVxCorrDelayMeasure(k)
% function LVxCorrDelayMeasure(k)
% Demonstrates the principle of identifying the time
% delay between two signals which are correlated but offset in time.
% A certain amount of white noise (amplitude set by the value of
% k is mixed into the signal.
% Typical Test call:
% LVxCorrDelayMeasure(0.1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

t = 0:1/32:1-1/32;
ysin = sin(4*pi*t);
y1(1,1:256)= zeros(1,256); 
y1(1,6:21)= ysin(1,1:16); % directly from the source
y2(1,1:512)= zeros(1,512); 
y2(1,41:56)= ysin(1,1:16);       % directly from the source
y2(1,81:96)= 0.7*ysin(1,1:16);   % a first echo
y2(1,111:126)= 0.4*ysin(1,1:16); % a second echo

thecorr(1,1:410) = zeros(1,410);

aa = 0:1:150;

ynoise = k*randn(1,575);
y1 = y1 + ynoise(1,1:256);
y2 = y2(1,1:512) + ynoise(1,64:575);

maxplot = 20;

figure(29);
clf

subplot(311)
stem(y1(1,1:256),'b.')
ylabel(['Amplitude'])
xlabel(['(a)  Samples from First Microphone; X-Axis = Sample Number'])
axis([0 256 -2 2])

for Lag = 0:1:150
   
subplot(312)
xvec = Lag:1:Lag+255;
stem(xvec,y2(1,1+Lag:256+Lag),'b.')
ylabel(['Amplitude'])
xlabel(['(b)  Samples from 2nd (further) Mic, shifted to Left ',...
      num2str(Lag),' Samples; X-Axis = Samp No.'])

thecorr(Lag+1) = sum(y1(1,1:256).*y2(1,1+Lag:256+Lag));
axis([Lag Lag+256 -2 2])

subplot(313)
stem(aa,thecorr(1,1:length(aa)),'b.')
ylabel(['Amplitude'])
xlabel(['(c)  The Correlation Sequence; X-Axis = Lag Number'])
axis([0 150 -12 12])

pause(0.02)
end

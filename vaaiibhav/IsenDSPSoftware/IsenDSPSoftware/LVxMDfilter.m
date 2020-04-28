function LVxMDfilter(k,lenMDfilt,plotlim)
% Creates a test signal sampled at 10 kHz consisting of, over a duration of 
% one second, fifty equally spaced bipolar pulses each consisting of a one 
% millisecond long positive pulse followed immediately by a one millisecond
% long negative pulse, the entire one second test signal consisting of the fifty 
% bipolar pulses plus white noise having standard deviation equal to k.
% A matched FIR filter of length lenMDfilt is used to improve the 
% signal-to-noise ratio, i.e., to emphasize the signal relative to the noise.
% plotlim is the number of samples of tstSig and the two filtered sequences to plot
% Test calls:
% LVxMDfilter(0.1,20,2000)
% LVxMDfilter(0.25,20,2000)
% LVxMDfilter(0.5,20,2000)
% LVxMDfilter(1,20,2000)
% LVxMDfilter(2,20,2000)
% LVxMDfilter(0.1,6,2000)
% LVxMDfilter(0.25,6,2000)
% LVxMDfilter(0.5,6,2000)
% LVxMDfilter(1,6,2000)
% LVxMDfilter(2,6,2000)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
col = [ones(1,10),-ones(1,10),zeros(1,180)]';
tstSig = col*ones(1,50);
tstSig = tstSig(:)';
tstSig = tstSig + k*randn(1,length(tstSig));

if lenMDfilt < 2
lenMDfilt = 2;
end
lenFiltov2 = fix(lenMDfilt/2);
filtLen = 2*lenFiltov2;
imp = fliplr([ones(1,lenFiltov2),-ones(1,lenFiltov2)]);

% first method, 
ans1 = conv(imp,tstSig);

% second method

ans2 = zeros(1,length(tstSig)+length(imp));
ans2(1,2*lenFiltov2) = sum(tstSig(1,1:lenFiltov2)) - sum(tstSig(1,lenFiltov2+1:filtLen));
for ctr = filtLen+1:1:length(tstSig)
ans2(1,ctr) = ans1(1,ctr-1) - tstSig(1,ctr) - tstSig(1,ctr-filtLen) + 2*tstSig(1,ctr-lenFiltov2);
end

plotlim = 2000;

figure(786)
clf

subplot(311)
plot(tstSig(1,1:plotlim))
ylabel('Amp')
xlabel('(a0 Sample')

subplot(312)
plot(ans1(1,1:plotlim))
ylabel('Amp')
xlabel('(b) Sample')

subplot(313)
plot(ans2(1,1:plotlim))
ylabel('Amp')
xlabel('(c) Sample')
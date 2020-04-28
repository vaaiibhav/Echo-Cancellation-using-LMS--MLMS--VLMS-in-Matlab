function LVxFIRFreqRespMultMeth(imp,lenEval)
% function LVxFIRFreqRespMultMeth([1,0,0,0,1],1024)
% Receives a real or complex impulse response and computes 
% and displays the frequency response using four methods, namely
% the DTFT, the z-transform, real chirp response, and complex
% chirp response. The magnitude of all four responses is plotted 
% on a single figure.
% imp is the impulse response to be evaluated
% lenEval is the number of frequency samples to compute
% All four frequency response methods test the frequency
% response from 0 to + 2*pi radians.
% Test calls:
% LVxFIRFreqRespMultMeth(ones(1,8),1024)
% LVxFIRFreqRespMultMeth(-ones(1,8),1024)
% LVxFIRFreqRespMultMeth([1,zeros(1,6),1],1024)
% LVxFIRFreqRespMultMeth([-1,zeros(1,6),-1],1024)
% LVxFIRFreqRespMultMeth( (ones(1,8) + j*[1,zeros(1,6),1]),1024)
% LVxFIRFreqRespMultMeth((ones(1,8) - j*[1,zeros(1,6),1]),1024)
% LVxFIRFreqRespMultMeth(exp(j*2*pi*[0:1:7]*3/8),1024)
% LVxFIRFreqRespMultMeth(exp(-j*2*pi*[0:1:7]*3/8),1024)
% LVxFIRFreqRespMultMeth([exp(-j*2*pi*[0:1:31]*5/32)+...
%  0.5*exp(j*2*pi*[0:1:31]*11/32)],1024)
% LVxFIRFreqRespMultMeth([cos(2*pi*[0:1:31]*5/32)+...
%  0.5*sin(2*pi*[0:1:31]*11/32)],1024)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

figure(988)

%if rem(lenEval,2)==0
%    loLim = -lenEval/2+1
%    hiLim = lenEval/2;
%else
%    loLim = -(lenEval-1)/2;
%    hiLim = (lenEval-1)/2;
%end

loLim = 0;
hiLim = lenEval;

subplot(221) % DTFT via long FFT
y = abs(fft(imp,lenEval)); 
%y = fftshift(y);
xvec = 2*[loLim:1:hiLim-1]/hiLim;
plot(xvec,y)
xlabel('(a) Normalized Frequency')
ylabel('Magnitude')

subplot(222) 
B = 0:1:length(imp)-1;
pzMat = (2*pi*[loLim:1:hiLim-1]/hiLim)'*(-B); 
ZZMat = exp(j*pzMat); 
zXform = ZZMat*conj(imp'); 
plot(xvec,abs(zXform));
xlabel('(b) Normalized Frequency')
ylabel('Magnitude')

subplot(223)
xvec = 2*[loLim:1:hiLim-1]/hiLim;
t = [0:1/(lenEval-1):1];
y = abs(filter(imp,1,chirp(t,loLim,1,hiLim)));
plot(xvec,y);
xlabel('(c) Normalized Frequency')
ylabel('Magnitude')

subplot(224)
ts = chirp(t,loLim,1,hiLim) - j*chirp(t,loLim,1,hiLim,'linear',90);
plot(xvec,abs(filter(imp,1,ts)))
xlabel('(d) Normalized Frequency')
ylabel('Magnitude')
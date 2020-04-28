% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% this script is not meant to be called; it contains clean copies of m-code
% from the text. To run any of the example m-code segments below, copy the
% code, paste it into the Command line, and press Return.
return
% Page 51------------------------------------------------------------------
LVxLMSANCNorm([0, 0, 1, 0, -0.5, 0.6, 0, 0, -1.2, 0],2,2,3,[1,2,4,8])
LVxLMSANCNorm([0,0,1,0,-0.5,0.6,0,0,-1.2,0],0,2,3,[1,2,4,8])
% Page 53------------------------------------------------------------------
LVxModelPlant([1,0,0.81],1,100,3,3,0.5,0,0,1000)
% Page 56------------------------------------------------------------------
LVxLMSAdptFiltEchoShort(0.01,0.2,1,1,50)

global EchoErr
sound(EchoErr,8000)
% Page 58------------------------------------------------------------------
LVxLMSAdaptFiltEcho(0.02,0.2,0,0,0,50,0.02)
LVxLMSAdptFiltEchoShort(0.02,0.2,0,0,50)
LVxLMSAdaptFiltEcho(0.02,0.2,0,0,1,50,0.02)
% Page 60------------------------------------------------------------------
LVxLMSAdaptFiltEcho(0.02,0.2,1,1,[],[],[])
% Page 63------------------------------------------------------------------
LVxLMSAdaptFiltDeCorr(0.2,1,100,10,2)
MLxLMSAdaptFiltDeCorr(0.2,1,100,10,2)
% Page 65------------------------------------------------------------------
LVxLMSAdaptFiltDeCorr(0.2,1,100,10,0)
% Page 69------------------------------------------------------------------
LVxLMSInterferCancel(0.02,0.3,1,1,0,50,0.03,1/16,1/16,6,24,6,24,18)
LVxLMSInterferCancShort(0.02,0.3,1,1,1/16,1/16,6,24,6,24,18,30)
LVxLMSInterferCancel(0.02,0.3,1,1,0,50,0.03,1,1,1,6,1,6,18)
% Page 72------------------------------------------------------------------
LVxChannEqualiz([1,0.7],2.15,2,17)
% Page 74------------------------------------------------------------------
LVxChannEqualiz([1,0,1],2.4,2,10)

y = conv([1,0,1],[0,1,0,-1,0,1,0,-1,0,1])
% y = [1,0,0,0,0,0,0,0,0,0,1]
fr = abs(fft(y,1024)); plot(fr(1,1:513))
% Page 76------------------------------------------------------------------
LVxReverbDereverb(1450,0.7,0.05,400) 

LVxRevDerevShort(325,0.7,0.04,200,30)
% Page 78------------------------------------------------------------------
% (use the following code to supply a reverberated audio signal as 
% ReverbSnd to run the code from page 78 below
Delay = 600; DecayRate = 0.8; PeakSep = 300;
[ReverbSnd,Fs,bits] = wavread('drwatsonSR4K.wav');
ReverbSnd = filter(1,[1 zeros(1,Delay-1) -DecayRate],ReverbSnd);

b = xcorr(ReverbSnd);  % this code segment from page 78
[locabsPeaks, absvalPeaks] = LVfindPeaks(abs(b),2,PeakSep);
absestDelay = abs(locabsPeaks(1) - locabsPeaks(2))



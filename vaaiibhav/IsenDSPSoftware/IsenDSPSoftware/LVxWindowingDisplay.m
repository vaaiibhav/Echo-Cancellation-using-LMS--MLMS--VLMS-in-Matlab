function LVxWindowingDisplay(NoiseAmp,freq1,freq2)
% LVxWindowingDisplay(0,66.6,70.5)
% LVxWindowingDisplay(2,66.6,70.5)
% LVxWindowingDisplay(0,60,70)
% LVxWindowingDisplay(2,60,70)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
figure(16167)
clf

SR = 256;
t = 0:1/SR:1-1/SR;
Amp1 = 1;
Amp2 = 1;
for WinType=1:1:4  
   if WinType==1
      w = rectwin(SR);
      stW = 'Rectangular';
      strMarker = '-';
   elseif WinType==2
      w = kaiser(SR,5);
      stW = 'Kaiser';
      strMarker = ':';
   elseif WinType==3
      w = blackman(SR);
      stW = 'Blackman';
      strMarker = '-.';
   elseif WinType==4
      w = hamming(SR);
      stW = 'Hamming';
      strMarker = '--';
   elseif WinType==5
      w = rectwin(SR);
      stW = 'Rectangular';
   end
   
SignalPart = Amp1*sin(2*pi*t*freq1) + Amp2*cos(2*pi*t*freq2);
netAbsFFT = zeros(1,SR);

for ctr = 1:1:50
y = SignalPart + NoiseAmp*randn(1,length(t));
ywin = (w').*y;
absywinfft = abs(fft(ywin));
netAbsFFT = netAbsFFT + absywinfft;
end

maxabsywinfft = max(netAbsFFT);
absywinfft = (1/maxabsywinfft)*netAbsFFT;
hold on
lowFreqDisp = max([0  freq1-5]);
hiFreqDisp = min([freq2+5  SR/2]);
plot(floor(lowFreqDisp):1:ceil(hiFreqDisp),20*log10(absywinfft(1,floor(lowFreqDisp)+1:ceil(hiFreqDisp)+1)),['b',strMarker]);
line([freq1 freq1],[-20 0]);
line([freq2 freq2],[-20 0]);

ylabel(['Magnitude, dB'])
xlabel(['FFT Bin Number; Rect=Solid; Kaiser(5)=Dotted; Blackman=DashDot; Hamming=Dashed'])
axis([lowFreqDisp hiFreqDisp  -30    0])
end
function LVSynthPlotSquareSawtooth(WaveType,FundFreq,SR)
% LVSynthPlotSquareSawtooth(1,10,1000)
%
%WaveType: 1 for sawtooth, 2 for squarewave, 3 for triangle
%FundFreq is a real number
%SR, sample rate, is a real positive integer
%The number of harmonics computed will be adjusted to be the maximum for
%the values of FundFreq and SR.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Inc = WaveType;
t = 0:1/SR:1-1/SR;
Waveform = zeros(1,SR);
lim = floor(0.5*SR);
x = 1:1:SR;

figure(900)
clf

subplot(2,1,1)
  xlabel(['(a)  Component Waves (X-axis = time, sec): Press Any Key to Continue']);
   ylabel(['Amplitude']);

   if Inc==1|Inc==2
      ThisInc=Inc;
   elseif Inc==3
      ThisInc=2;
   else
      Comment='Waveform type must be entered as 1, 2, or 3; ending'
      return
   end
   
for n = 1:ThisInc:fix(abs(lim/FundFreq))
   
   phi = 0;
 %  phi = rand(1);
   if Inc==1|Inc==2
      LatestHarmonic = (1/n)*sin(2*pi*n*t*FundFreq + phi );
   elseif Inc==3
      LatestHarmonic = (1/n^2)*cos(2*pi*n*t*FundFreq + phi );  % Triangle
   end
      
   Waveform = Waveform + LatestHarmonic;
   
   subplot(2,1,1)
   hold on
   plot((x-1)/SR,LatestHarmonic)
   % stem(x-1,LatestHarmonic,'b.')
   axis([0 (SR-1)/SR -1.2 1.2])
   hold off
   
   subplot(2,1,2)
   plot((x-1)/SR,Waveform,'b')
   xlabel(['(b)  Summation Of All Components up to Harmonic ',num2str(n)]);
   ylabel(['Amplitude']);

   plmax = 1.2*max(abs(Waveform));
   axis([0 (SR-1)/SR -plmax plmax])
   pause
end
subplot(2,1,2)
hold on
plot((x-1)/SR,Waveform,'k')


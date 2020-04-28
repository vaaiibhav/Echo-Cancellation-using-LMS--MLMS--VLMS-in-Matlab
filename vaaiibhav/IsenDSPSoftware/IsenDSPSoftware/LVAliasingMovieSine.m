function LVAliasingMovieSine
%
%function LVAliasingMovieSine
%
%Opens up a GUI that automatically resamples eight cycles of a sinusoid at
%sample rates starting with 240 and decreasing to 8. When the sample rate
%gets near the Nyquist limit, the user is notified to press any key to
%continue computation.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

Inc = 200;

Freq = 8;

RefWaveform = sin(2*pi*Freq*(0:1/240:1));

SR = [240,140,80,40,30,24,20,17.5,17,16.5,16.1,16.01,16,15.99,15.9,15.5,15,14.5,14,13.5,13,12.5,12,11.5,11,10.5,10,9.8,9.6,9.4,9.2,9.0,8.9,8.8,8.7,8.6,8.5,8.4,8.3,8.2,8.1,8.08,8.06,8.04,8.02,8.001 8.0];

figure(51)

for Index = 1:1:length(SR);

Waveform = []; 

t = 0:1/SR(Index):1;

Waveform = sin(2*pi*Freq*t);

wavchar = 'Sine';

subplot(2,1,1)
plot((0:1:240)/240,RefWaveform,'k--')
hold on
stem(t,Waveform,'b.')
plot(t,Waveform,'r')
hold off
ylabel(['Amplitude'])
%xlabel(['(a)  A Sine Wave, Freq = 8 Hz (SR = ',num2str(SR(Index)),' Hz); X-Axis = Time, sec'],'fontsize',12)
 xlabel(['(a)  Time, seconds; SR = ',num2str(SR(Index)),' Hz'])

axis([0 inf -1.1 1.1])

subplot(2,1,2)

vbn = Waveform;
xxx = abs(fft(vbn));
 
stem(0:1:ceil(0.5*length(xxx))-1,xxx(1,1:ceil(0.5*length(xxx))),'ko');
ylabel(['Magnitude'])
xlabel(['(b)  Frequency Content of Sequence in Hz'])
axis([0 30 0 (5 + max(xxx))])

if SR(Index) < 24
 subplot(2,1,2)
   text(15,fix((5 + max(xxx))/2),'Press Any Key to Continue','color',[1 0 0]);
   pause
else
   pause(0.5)
end

end


subplot(2,1,2)
xlabel('Demo Finished!');




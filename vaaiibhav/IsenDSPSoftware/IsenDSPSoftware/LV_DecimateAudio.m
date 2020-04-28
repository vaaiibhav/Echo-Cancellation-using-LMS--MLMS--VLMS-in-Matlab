function LV_DecimateAudio(DecimateFac)
% global decsig origsig thisSR
% function LV_DecimateAudio(DecimateFac)
%
% DecimatFac is the factor by which to decimate
% a mixture of three sinusoids having original
% frequencies of 200, 400, and 600 Hz and
% may have a value of 2 or 4
%
% A typical call would be 
%
% LV_DecimateAudio(2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(DecimateFac==2|DecimateFac==4)
   error('DecimateFac must be 2 or 4')
end

DisplayLim = floor(900/DecimateFac);

SR = 44100;

t = 0:1/SR:1-1/SR;

origsig = zeros(1,SR);

for freq = 200:200:600
   origsig = origsig + sin(2*pi*t*freq);
end

origsig = (0.99/max(abs(origsig)))*origsig;

y = abs(fft(origsig,44100));

sound(origsig,44100)

cutoff = 1/DecimateFac;
a=1;
b = fir1(22,cutoff);

decsig = conv(b,origsig);

ind = 1:DecimateFac:length(decsig);
decsig = decsig(1,ind);

decsig = (0.99/max(abs(decsig)))*decsig;
start = floor(22/(DecimateFac+1)-1);
decsig = decsig(1,start:length(decsig));

lendecsig = length(decsig);
ftdecsig = abs(fft(decsig,44100));

pause(length(origsig)/44100)
sound(decsig,44100)

figure(2301)
clf

subplot(211)
stem(origsig(1,1:DisplayLim))
xlabel(['(a)  The Original Sequence (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0, DisplayLim, -inf, inf])

subplot(212)
plot(y(1,1:22050))
xlabel(['(b)  Frequency, Original Sequence'])
ylabel(['Magnitude'])
axis([0, 3000, -inf, inf])

figure(2302)
clf

subplot(211)
stem(decsig(1,1:DisplayLim))
xlabel(['(a)  Decimated Sequence (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0, DisplayLim, -inf, inf])

subplot(212)
plot(ftdecsig(1,1:22050))
xlabel(['(b)  Frequency, Decimated Seq'])
ylabel(['Magnitude'])
axis([0 3000 -inf inf])

figure(2303)
clf

subplot(2,1,1)
stem(origsig(1,1:DisplayLim))
xlabel(['(a)  The Original Sequence (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0 DisplayLim -inf inf])

subplot(2,1,2)
stem(decsig(1,1:DisplayLim))
xlabel(['(b)  The Decimated Sequence (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0 DisplayLim -inf inf])

figure(2304)
clf

subplot(2,1,1)
plot(y(1,1:22050))
xlabel(['(a)  Frequency, Original Sequence'])
ylabel(['Magnitude'])
axis([0 3000 -inf inf])

subplot(2,1,2)
plot(ftdecsig(1,1:22050))

xlabel(['(b)  Frequency, Decimated Sequence'])
ylabel(['Magnitude'])
axis([0 3000 -inf inf])


pause(length(decsig)/44100)
sound(decsig,44100/DecimateFac)




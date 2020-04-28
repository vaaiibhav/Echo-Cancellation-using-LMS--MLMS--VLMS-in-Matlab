function LV_InterpAudio(InterpFac)
% function LV_InterpAudio(InterpFac)
%
% InterpFac is an integer which is a ratio of the number of samples
% the test signal will have after interpolation to the original number;
% InterpFac may have values of 2 or 4.
%
% A typical call:
%
% LV_InterpAudio(4)
% LV_InterpAudio(2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

if ~(InterpFac==2|InterpFac==4)
   error('InterpFac must be either 2 or 4')
end

SR = 11025;

FinalSR = InterpFac*SR;
DisplayLim = 225;

t = 0:1/SR:1-1/SR;

OrigSig = zeros(1,SR);

for freq = 400:400:1200
   OrigSig = OrigSig + cos(2*pi*t*freq + pi/8);
end

aaa = max(abs(OrigSig));

OrigSig = (0.99/aaa)*OrigSig;

lenOrigSig = length(OrigSig);
sound(OrigSig,SR)

y = abs(fft(OrigSig,2^17));

ind = InterpFac:InterpFac:InterpFac*lenOrigSig;
UpSampSig = zeros(1,InterpFac*lenOrigSig);
UpSampSig(1,ind) = OrigSig(1,ind/InterpFac);

cutoff = 1/InterpFac;

a=1;

b = fir1(22,cutoff);

XUpSampSig = conv(b,UpSampSig);
UpSampSig = XUpSampSig(1,33:length(XUpSampSig));
XUpSampSig = [];

bbb = max(abs(UpSampSig));

UpSampSig = (0.99/bbb)*UpSampSig;
lenUpSampSig = length(UpSampSig);

ftUpSampSig = abs(fft(UpSampSig,2^17));

pause(lenOrigSig/SR)
sound(UpSampSig,SR)

figure(2001)
clf

subplot(211)
stem(OrigSig(1,1:DisplayLim));
xlabel(['(a)  The Original Sequence  (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0 DisplayLim -inf inf])

subplot(212)
plot((22050/4/2^16)*(1:1:2^16),y(1,1:2^16))
xlabel(['(b)  Frequency, Original Sequence'])
ylabel(['Magnitude'])
axis([0 3000 -inf inf])

figure(2002)
clf

subplot(211)
stem(UpSampSig(1,1:DisplayLim))
xlabel(['(a)  Upsampled Sequence  (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0 DisplayLim -inf inf])

subplot(212)
plot((22050/4/2^16)*(1:1:2^16),ftUpSampSig(1,1:2^16))
xlabel(['(b)  Frequency, Upsampled Seq'])
ylabel(['Magnitude'])
axis([0 3000 -inf inf])

figure(2003)
clf

subplot(2,1,1)
stem(OrigSig(1,1:DisplayLim));
xlabel(['(a)  Original Sequence  (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0 DisplayLim -inf inf])

subplot(2,1,2)
stem(UpSampSig(1,1:DisplayLim));
xlabel(['(b)  Upsampled Sequence; InterpFac = ',num2str(InterpFac),' (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0 DisplayLim -inf inf])

figure(2004)
clf
subplot(2,1,1)
plot((22050/4/2^16)*(1:1:2^16),y(1,1:2^16))
xlabel(['(a)  Frequency, Original Sequence'])
ylabel(['Magnitude'])
axis([0 3000 -inf inf])

subplot(2,1,2)
plot( (22050/4/2^16)*(1:1:2^16),ftUpSampSig(1,1:2^16))
xlabel(['(b) Frequency, Upsampled Sequence'])
ylabel(['Magnitude'])
axis([0 3000 -inf inf])

pause(InterpFac*lenOrigSig/SR)
sound(UpSampSig,FinalSR)

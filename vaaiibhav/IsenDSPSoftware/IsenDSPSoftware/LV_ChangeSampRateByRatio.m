function LV_ChangeSampRateByRatio(OrigSR,InterpFac,DecimateFac)
% function LV_ChangeSampRateByRatio(OrigSR,InterpFac,DecimateFac)
%
% OrigSR is the sample rate to use to create an original audio signal
% If OrigSR is any of 8000, 11025,22050,
% or 44100 Hz., the script will play the original audio signal
% If FinalSR as computed by the script from the input arguments is within 1% of any of 
% 8000, 11025,22050, or 44100 Hz., the script will play the resampled audio
% signal at the nearest standard sample rate
%
% Typical calls:
%
% LV_ChangeSampRateByRatio(22050,8,22)
% LV_ChangeSampRateByRatio(44100,8,44)
% LV_ChangeSampRateByRatio(8000,22,8)
% LV_ChangeSampRateByRatio(8000,44,8)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

SR = OrigSR;

FinalSR = (InterpFac/DecimateFac)*OrigSR
% DisplayLim = min([fix(200*InterpFac/DecimateFac) 500]);

DisplayLim = 225;

t = 0:1/SR:1-1/SR;

OrigSig = 0;

for freq = 400:400:1200
   OrigSig = OrigSig + sin(2*pi*t*freq + pi/8);
end

aaa = max(abs(OrigSig));

OrigSig = (0.99/aaa)*OrigSig;

lenOrigSig = length(OrigSig)
playOrig = 0;
if OrigSR==8000|OrigSR==11025|OrigSR==22050|OrigSR==44100
sound(OrigSig,OrigSR)
playOrig=1;
end

y = abs(fft(OrigSig,OrigSR));

%==========================================================
ind = 1:InterpFac:InterpFac*lenOrigSig;
UpSampSig = zeros(1,InterpFac*lenOrigSig);
UpSampSig(1,ind) = OrigSig(1,1:1:lenOrigSig);
cutoff = 0.999/max([InterpFac DecimateFac])
%==========================================================

a=1;

filtlen = 40 + 3*round(1/cutoff);
b = fir1(filtlen,cutoff);

UpSampSig = filter(b,a,UpSampSig);

bbb = max(abs(UpSampSig));

UpSampSig = (0.99/bbb)*UpSampSig;
lenUpSampSig = length(UpSampSig);

DecIndices = 1+fix(filtlen/2):DecimateFac:lenUpSampSig; % rough compensation for filter delay

ResampledSig = UpSampSig(1,DecIndices);

ftResampledSig = abs(fft(ResampledSig,OrigSR));
lenResampSig = length(ResampledSig)

UpSampSig = [];

if playOrig==1
    pause(2)
end

Pfrs = [8000,11025,22050,44100];
frPCDiffs = (Pfrs-FinalSR)/FinalSR;
bt = find(abs(frPCDiffs) < 0.01);
if ~(isempty(bt))
PlaySR = Pfrs(bt(1));
sound(ResampledSig,PlaySR)
end


figure(2201)
clf

subplot(211)
stem(OrigSig(1,1:DisplayLim));
xlabel(['(a)  Original Sequence (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0, DisplayLim, 1.1*min(OrigSig), 1.1*max(OrigSig)])

subplot(212)

plot(y(1,1:fix(length(y)/2)))
xlabel(['(b)  Frequency, Original Sequence'])
ylabel(['Magnitude'])
axis([0 OrigSR/2 -inf inf])

figure(2202)
clf

subplot(211)
stem(ResampledSig(1,1:DisplayLim));
xlabel(['(a)  Resampled Sequence (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0, DisplayLim, 1.1*min(ResampledSig), 1.1*max(ResampledSig)])

subplot(212)
plot(ftResampledSig(1,1:fix(length(ftResampledSig)/2)))
xlabel(['(b)  Frequency, Resampled Seq'])
ylabel(['Magnitude'])
axis([0 OrigSR/2 -inf inf])

figure(2203)
clf

subplot(2,1,1)
stem(OrigSig(1,1:DisplayLim));
xlabel(['(a)  Original Sequence  (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0, DisplayLim, 1.1*min(OrigSig), 1.1*max(OrigSig)])

subplot(2,1,2)
stem(ResampledSig(1,1:DisplayLim));
xlabel(['(b)  Resampled Sequence; InterpFac = ',num2str(InterpFac),' & DecFac = ',num2str(DecimateFac),' (1st ',num2str(DisplayLim),' Samples)'])
ylabel(['Amplitude'])
axis([0, DisplayLim, 1.1*min(ResampledSig), 1.1*max(ResampledSig)])


figure(2204)
clf

subplot(2,1,1)
plot( y(1,1:fix(length(y)/2)))
xlabel(['(a)  Frequency, Original Sequence'])
ylabel(['Magnitude'])
axis([0 OrigSR/2 -inf inf])

subplot(2,1,2)
plot(ftResampledSig(1,1:fix(length(ftResampledSig)/2)))
xlabel(['(b)  Frequency, Resampled Seq'])
ylabel(['Magnitude'])
axis([0 OrigSR/2 -inf inf])




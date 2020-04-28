function ML_ChangeSampRateByRatioPOC(InterpFac,DecimateFac)
global ResampledSig2201 TheOrigSig2201 FinalSR2201 OrigSR2201
% function ML_ChangeSampRateByRatioPOC(InterpFac,DecimateFac)
%
% InterpFac is an integer that is a ratio of the number of samples
% the test signal will have after interpolation to the original number;
%
% InterpFac=2, for example, means there will be twice as many samples after
% interpolation as before.
%
% If DecimateFac = n, then the decimation process is performed by
% taking every nth sample of the interpolated sequence to send to the (final)
% output sequence.
%
% Appropriate filtering is performed (see text)
%
% A typical call:
%
% ML_ChangeSampRateByRatioPOC(3,4)
% ML_ChangeSampRateByRatioPOC(2,7)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
scnsize = get(0,'ScreenSize');

pos3 = [0.0025*scnsize(3), 0.05*scnsize(4), 0.995*scnsize(3), 0.84*scnsize(4)];
pos1 = pos3;
pos2 = pos3;
pos2a = pos3;

JFac = 1;
OrigSR2201 = JFac*44100;
SR = OrigSR2201;

FinalSR2201 = (InterpFac/DecimateFac)*OrigSR2201;
DisplayLim = min([fix(200*InterpFac/DecimateFac) 500]);

t = 0:1/SR:1-1/SR;

TheOrigSig2201 = zeros(1,SR);

for freq = 3000:1000:5000
   TheOrigSig2201 = TheOrigSig2201 + sin(2*pi*t*freq + pi/8);
end

aaa = max(abs(TheOrigSig2201));

TheOrigSig2201 = (0.99/aaa)*TheOrigSig2201;

lenOrigSig2201 = length(TheOrigSig2201);

y = abs(fft(TheOrigSig2201,2^17));

%==========================================================
ind = 1:InterpFac:InterpFac*lenOrigSig2201;
UpSampSig = zeros(1,InterpFac*lenOrigSig2201);
UpSampSig(1,ind) = TheOrigSig2201(1,1:1:lenOrigSig2201);
cutoff = 0.999/max([InterpFac DecimateFac]);
%==========================================================

a=1;
filtlen = 40;
b = fir1(filtlen,cutoff);

UpSampSig = filter(b,a,UpSampSig);

bbb = max(abs(UpSampSig));

UpSampSig = (0.99/bbb)*UpSampSig;
lenUpSampSig = length(UpSampSig);

DecIndices = 1+fix(filtlen/2):DecimateFac:lenUpSampSig; % rough compensation for filter delay

ResampledSig2201 = UpSampSig(1,DecIndices);

ftResampledSig = abs(fft(ResampledSig2201,2^17));

figure(2201)
clf
set(2201,'Color',[1 1 1],'Position',pos2,'numbertitle','off','name','Original Sequence and Spectrum')

haxSP650Upper = axes('position',[0.15 0.6 0.8 0.3]);
haxSP650Lower = axes('position',[0.15 0.15 0.8 0.3]);

uicontrol('Style','push','Callback','PlayOrigSRSeq2201','String','PlayOrigSequence',...
   'Position',[5 4 100 19],'backgroundcolor',[1 1 1]);

subplot(haxSP650Upper)
plot(TheOrigSig2201(1,1:DisplayLim))
xlabel(['(a)  The Original Sequence (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(haxSP650Lower)
plot(y(1,1:2^16))
xlabel('(b)  Spectrum of the Original Sequence (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 2^16 -inf inf])

figure(2202)
clf
set(2202,'Color',[1 1 1],'Position',pos2a,'numbertitle','off','name','Resampled Sequence and Spectrum')

haxSP650Upper = axes('position',[0.15 0.6 0.8 0.3]);
haxSP650Lower = axes('position',[0.15 0.15 0.8 0.3]);

subplot(haxSP650Upper)
plot(ResampledSig2201(1,1:DisplayLim))
xlabel(['(a)  The Resampled Sequence (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(haxSP650Lower)
plot(ftResampledSig(1,1:2^16))
xlabel('(b)  Spectrum of Resampled Seq (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 2^16 -inf inf])

uicontrol('Style','push','Callback','PlayResampledSig2201','String','Play at New Samp Rate',...
   'Position',[225 4 150 19],'backgroundcolor',[1 1 1]);

uicontrol('Style','push','Callback','PlayResampSigOR2201','String','Play Resamp Seq at Orig Samp Rate',...
   'Position',[10 4 200 19],'backgroundcolor',[1 1 1]);

figure(2203)
clf
set(2203, 'position',pos1,'color',[1 1 1],'numbertitle','off','name','Original and Resampled Sequences')

subplot(2,1,1)
plot(TheOrigSig2201(1,1:DisplayLim))
set(gca,'fontsize',18)
xlabel('(a)  Sample Number','fontsize',18)
%xlabel(['(a)  The Original Sequence  (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',18)
axis([0 DisplayLim -inf inf])

subplot(2,1,2)
plot(ResampledSig2201(1,1:DisplayLim))
set(gca,'fontsize',18)
xlabel('(b)  Sample Number','fontsize',18)
%xlabel(['(b)  The Resampled Sequence; InterpFac = ',num2str(InterpFac),' & DecFac = ',num2str(DecimateFac),' (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',18)
axis([0 DisplayLim -inf inf])

%wavwrite(TheOrigSig,OrigSR,8,['DemoChangeSampRateByRatioOrigSeqInterp',num2str(InterpFac),'Dec',num2str(DecimateFac)])

%wavwrite(ResampledSig,FinalSR,8,['DemoChangeSampRateByRatioResampSeqInterp',num2str(InterpFac),'Dec',num2str(DecimateFac)])

%wavwrite(ResampledSig,OrigSR,8,['DemoChangeSampRateByRatioResampSeqOSRInterp',num2str(InterpFac),'Dec',num2str(DecimateFac)])

figure(2204)
clf
set(2204, 'position',pos3,'color',[1 1 1],'numbertitle','off','name','Spectra of Original and Resampled Sequences')

subplot(2,1,1)
plot(y(1,1:2^16))
xlabel('(a)  Spectrum of the Original Sequence (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 2^16 -inf inf])

subplot(2,1,2)
plot(ftResampledSig(1,1:2^16))
xlabel('(b)  Spectrum of Resampled Sequence (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 2^16 -inf inf])




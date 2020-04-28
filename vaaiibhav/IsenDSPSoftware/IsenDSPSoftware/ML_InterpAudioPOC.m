function ML_InterpAudioPOC(InterpFac)
global ResampledSig2001 TheOrigSig2001 FinalSR2001 OrigSR2001
%function ML_InterpAudioPOC(InterpFac)
%
%InterpFac is an integer which is a ratio of the number of samples
%the test signal will have after interpolation to the original number;
%InterpFac=2, for example, means there will be twice as many samples after
%interpolation as before. Appropriate filtering is performed (see text)
%
%A typical call:
%
%ML_InterpAudioPOC(3)
%ML_InterpAudioPOC(2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if InterpFac==1
   InterpFac = 2
end

DecimateFac = 1;

scnsize = get(0,'ScreenSize');
pos1 = [0.0025*scnsize(3), 0.05*scnsize(4), 0.995*scnsize(3), 0.84*scnsize(4)];
pos2 = [0.0025*scnsize(3), 0.05*scnsize(4), 0.995*scnsize(3), 0.88*scnsize(4)];
%pos2a = [0.005*scnsize(3), 0.05*scnsize(4), 0.99*scnsize(3), 0.8*scnsize(4)];
%pos3 = [0.005*scnsize(3), 0.05*scnsize(4), 0.99*scnsize(3), 0.85*scnsize(4)];

JFac = 1;
OrigSR2001 = JFac*11025;
SR = OrigSR2001;

FinalSR2001 = (InterpFac/DecimateFac)*OrigSR2001;
DisplayLim = min([fix(200*InterpFac/DecimateFac) 500]);

t = 0:1/SR:1-1/SR;

TheOrigSig2001 = zeros(1,SR);

for freq = 3000:1000:5000
   TheOrigSig2001 = TheOrigSig2001 + sin(2*pi*t*freq + pi/8);
end

aaa = max(abs(TheOrigSig2001));

TheOrigSig2001 = (0.99/aaa)*TheOrigSig2001;

lenOrigSig = length(TheOrigSig2001);

y = abs(fft(TheOrigSig2001,2^17));

ind = InterpFac:InterpFac:InterpFac*lenOrigSig;
UpSampSig = zeros(1,InterpFac*lenOrigSig);
UpSampSig(1,ind) = TheOrigSig2001(1,ind/InterpFac);

cutoff = 1/max([InterpFac DecimateFac]);

a=1;
b = fir1(62,cutoff);

XUpSampSig = conv(b,UpSampSig);
UpSampSig = XUpSampSig(1,33:length(XUpSampSig));
XUpSampSig = [];

bbb = max(abs(UpSampSig));

UpSampSig = (0.99/bbb)*UpSampSig;
lenUpSampSig = length(UpSampSig);

ftUpSampSig = abs(fft(UpSampSig,2^17));

ResampledSig2001 = UpSampSig;

figure(2001)
clf
set(2001,'Color',[1 1 1],'Position',pos2,'numbertitle','off','name','The Original Sequence and Its Spectrum')

haxSP650Upper = axes('position',[0.15 0.6 0.8 0.3]);
haxSP650Lower = axes('position',[0.15 0.15 0.8 0.3]);

PlOrigHF = uicontrol('Style','push','Callback','PlayOrigSRSeq2001','String','PlayOrigSequence',...
   'Position',[5 4 100 19],'backgroundcolor',[1 1 1]);

subplot(haxSP650Upper)
plot(TheOrigSig2001(1,1:DisplayLim))
xlabel(['(a)  The Original Sequence  (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
%set(gca,'fontsize',18)
%xlabel(['(a)  Sample Number'],'fontsize',12)
ylabel(['Amplitude'],'fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(haxSP650Lower)
plot((22050/2^16)*(1:1:2^16),y(1,1:2^16))
%set(gca,'fontsize',18)
%xlabel(['(b)  Sample Number'],'fontsize',12)
xlabel(['(b)  Spectrum of the Original Sequence (X-Axis = Bin Number)'],'fontsize',12)
ylabel(['Magnitude'],'fontsize',12)
axis([0 22050 -inf inf])

figure(2002)
clf
set(2002,'Color',[1 1 1],'Position',pos2,'numbertitle','off','name','The Upsampled (Interpolated) Sequence and Its Spectrum')

haxSP650Upper = axes('position',[0.15 0.6 0.8 0.3]);
haxSP650Lower = axes('position',[0.15 0.15 0.8 0.3]);

subplot(haxSP650Upper)
plot(UpSampSig(1,1:DisplayLim))
%set(gca,'fontsize',18)
%xlabel(['(a)  Sample Number'],'fontsize',18)
xlabel(['(a)  The Upsampled Sequence  (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel(['Amplitude'],'fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(haxSP650Lower)
plot((22050/2^16)*(1:1:2^16),ftUpSampSig(1,1:2^16))
%set(gca,'fontsize',18)
%xlabel(['(b)  Sample Number'],'fontsize',18)
xlabel(['(b)  Spectrum of Upsampled Seq (X-Axis = Bin Number)'],'fontsize',12)
ylabel(['Magnitude'],'fontsize',12)
axis([0 22050 -inf inf])

PlNewLF = uicontrol('Style','push','Callback','PlayResampledSig2001','String','Play at New Samp Rate',...
   'Position',[225 4 150 19],'backgroundcolor',[1 1 1]);

PlNewLFOrigRate = uicontrol('Style','push','Callback','PlayResampSigOR2001','String','Play Upsamp Seq at Orig Samp Rate',...
   'Position',[10 4 200 19],'backgroundcolor',[1 1 1]);

figure(2003)
clf
set(2003, 'position',pos1,'color',[1 1 1],'numbertitle','off','name','The Original and Upsampled (Interpolated) Sequences')

subplot(2,1,1)
plot(TheOrigSig2001(1,1:DisplayLim))
%set(gca,'fontsize',18)
%xlabel(['(a)  Sample Number'],'fontsize',18)
xlabel(['(a)  The Original Sequence  (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel(['Amplitude'],'fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(2,1,2)
plot(UpSampSig(1,1:DisplayLim))
%set(gca,'fontsize',18)
%xlabel(['(b)  Sample Number'],'fontsize',18)
xlabel(['(b)  The Upsampled Sequence; InterpFac = ',num2str(InterpFac),' (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel(['Amplitude'],'fontsize',12)
axis([0 DisplayLim -inf inf])

%wavwrite(TheOrigSig2001,OrigSR2001,8,['DemoInterpAudioInterpFac',num2str(InterpFac),'OrSeq'])

%wavwrite(UpSampSig,FinalSR2001,8,['DemoInterpAudioInterpFac',num2str(InterpFac),'FinSR'])

%wavwrite(UpSampSig,OrigSR2001,8,['DemoInterpAudioInterpFac',num2str(InterpFac),'OrigSR2001'])

figure(2004)
clf
set(2004, 'position',pos1,'color',[1 1 1],'numbertitle','off','name','Spectra of Original and Upsampled (Interpolated) Signals')

subplot(2,1,1)
plot((22050/2^16)*(1:1:2^16),y(1,1:2^16))
%set(gca,'fontsize',18)
%xlabel(['(a)  Bin Number'],'fontsize',18)
xlabel(['(a)  Spectrum of the Original Sequence (X-Axis = Bin Number)'],'fontsize',12)
ylabel(['Magnitude'],'fontsize',12)
axis([0 22050 -inf inf])

subplot(2,1,2)
plot( (22050/2^16)*(1:1:2^16),ftUpSampSig(1,1:2^16))
%set(gca,'fontsize',18)
%xlabel(['(b)  Bin Number'],'fontsize',18)
xlabel(['(b)  Spectrum of Upsampled Sequence (X-Axis = Bin Number)'],'fontsize',12)
ylabel(['Magnitude'],'fontsize',12)
axis([0 22050 -inf inf])



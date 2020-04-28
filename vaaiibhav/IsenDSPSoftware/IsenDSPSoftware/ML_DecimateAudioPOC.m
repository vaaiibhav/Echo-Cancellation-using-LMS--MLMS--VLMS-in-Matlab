function ML_DecimateAudioPOC(DecimateFac)
global decsig origsig thisSR
%function ML_DecimateAudioPOC(DecimateFac)
%
%DecimatFac is the factor by which to decimate
%a mixture of three sinusoids having original
%frequencies of 1000, 2000, and 3000 Hz
%
%A typical call would be 
%
%ML_DecimateAudioPOC(2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if DecimateFac==1
   DecimateFac = 2
end

DisplayLim = floor(400/DecimateFac);

scnsize = get(0,'ScreenSize');

pos1 = [0.0025*scnsize(3), 0.05*scnsize(4), 0.995*scnsize(3), 0.84*scnsize(4)];
pos2 = [0.0025*scnsize(3), 0.05*scnsize(4), 0.995*scnsize(3), 0.89*scnsize(4)];
pos2a = [0.0025*scnsize(3), 0.05*scnsize(4), 0.995*scnsize(3), 0.89*scnsize(4)];
pos3 = [0.0025*scnsize(3), 0.05*scnsize(4), 0.995*scnsize(3), 0.84*scnsize(4)];

JFac = 1;
SR = JFac*44100;
thisSR = SR;

t = 0:1/SR:1-1/SR;

origsig = zeros(1,SR);

for freq = 1000:1000:3000
   origsig = origsig + sin(2*pi*t*freq);
end

origsig = (0.99/max(abs(origsig)))*origsig;

y = abs(fft(origsig,44100));

cutoff = 1/DecimateFac;

b = fir1(22,cutoff);

decsig = conv(b,origsig);

ind = 1:DecimateFac:length(decsig);
decsig = decsig(1,ind);

decsig = (0.99/max(abs(decsig)))*decsig;
start = floor(22/(DecimateFac+1)-1);
decsig = decsig(1,start:length(decsig));

ftdecsig = abs(fft(decsig,44100));

figure(2301)
clf
set(2301,'Color',[1 1 1],'Position',pos2,'numbertitle','off','name','The Original Sequence and Its Spectrum')

haxSP650Upper = axes('position',[0.15 0.6 0.8 0.3]);
haxSP650Lower = axes('position',[0.15 0.15 0.8 0.3]);

uicontrol('Style','push','Callback','PlayHFSeq','String','PlayOrigSequence',...
   'Position',[5 4 100 19],'backgroundcolor',[1 1 1]);

subplot(haxSP650Upper)
plot(origsig(1,1:DisplayLim))
xlabel(['(a)  The Original Sequence (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(haxSP650Lower)
plot(y(1,1:22050))
xlabel('(b)  Spectrum of the Original Sequence (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 22050 -inf inf])

figure(2302)
clf
set(2302,'Color',[1 1 1],'Position',pos2a,'numbertitle','off','name','The Decimated Sequence and Its Spectrum')

haxSP650Upper = axes('position',[0.15 0.6 0.8 0.3]);
haxSP650Lower = axes('position',[0.15 0.15 0.8 0.3]);

subplot(haxSP650Upper)
plot(decsig(1,1:DisplayLim))
xlabel(['(a)  The Decimated Sequence (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(haxSP650Lower)
plot(ftdecsig(1,1:22050))
xlabel('(b)  Spectrum of Decimated Seq (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 22050 -inf inf])

uicontrol('Style','push','Callback','PlayLFSeq','String','PlayDecSeq',...
   'Position',[5 4 75 19],'backgroundcolor',[1 1 1]);

figure(2303)
clf
set(2303, 'position',pos1,'color',[1 1 1],'numbertitle','off','name','Original and Decimated Sequences')

subplot(2,1,1)
plot(origsig(1,1:DisplayLim))
xlabel(['(a)  The Original Sequence (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',12)
axis([0 DisplayLim -inf inf])

subplot(2,1,2)
plot(decsig(1,1:DisplayLim))
xlabel(['(b)  The Decimated Sequence (1st ',num2str(DisplayLim),' Samples)'],'fontsize',12)
ylabel('Amplitude','fontsize',12)
axis([0 DisplayLim -inf inf])

%wavwrite(origsig,thisSR,8,['DemoDecimateAudioOrigSeq'])

%wavwrite(decsig,thisSR,8,['DemoDecimateAudioDecFac',num2str(DecimateFac)])

figure(2304)
clf
set(2304, 'position',pos3,'color',[1 1 1],'numbertitle','off','name','Spectra of Original and Decimated Sequences')

subplot(2,1,1)
plot(y(1,1:22050))
xlabel('(a)  Spectrum of Original Sequence (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 22050 -inf inf])

subplot(2,1,2)
plot(ftdecsig(1,1:22050))

xlabel('(b)  Spectrum of Decimated Sequence (X-Axis = Bin Number)','fontsize',12)
ylabel('Magnitude','fontsize',12)
axis([0 22050 -inf inf])




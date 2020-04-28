function LV_InterpGeneric(Freq,InterpFac)
% function LV_InterpGeneric(4,2)
% This program inserts N-1 zeros between each
% sample of a set of samples and lowpass filters
% the result. The built-in MathScript "interp"
% function is used to do the same thing and is also 
% plotted. 
%
% Freq may be any positive integer 
% (16 being the maximum non-aliased frequency);
%
% InterpFac may be a positive integer such as 2,3,5,7, etc.
% A typical call:
% LV_InterpGeneric(8,2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

SR = 32;
t = 0:1/SR:1-1/SR;
y = sin(2*pi*t*Freq);
Newy(1,1:SR*InterpFac) = 0
n = 1:1:SR;
Newy(1,InterpFac*(n-1)+1)= y(1,n);

figure(33)

subplot(311)
stem(t*SR,y,'bo');
ampy = max(abs(y));
ylabel(['Amp'])
xlabel(['(a)  Sample, Original Sequence'])
axis([0 SR -1.2*ampy 1.2*ampy])

subplot(312)
stem(Newy,'bo');
ampNewy = max(abs(Newy));
ylabel(['Amp'])
xlabel(['(b)  Orig. Seq. w/ Zeros Between Samples'])
axis([0 InterpFac*SR -1.2*ampNewy 1.2*ampNewy])

cutoff = 1/InterpFac;
a=1;
b = fir1(22,cutoff);

% XOutputInterp = filter(b,a,Newy);
XOutputInterp = conv(b,Newy);

OutputInterp = [XOutputInterp(1,12:length(XOutputInterp)) ];
XOutputInterp = [];

plotlim = InterpFac*1.2*max(abs(OutputInterp));

subplot(313)
stem(InterpFac*OutputInterp,'bo');
ylabel(['Amp'])
xlabel(['(c)  The LP-Filtered, Upsampled Sequence'])
axis([0 InterpFac*SR -plotlim plotlim])

figure(34)

subplot(211)
stem(t*SR,y,'bo');
ylabel(['Amplitude'])
xlabel(['(a)  The Original Seq.'])
axis([0 SR -1.2*ampy 1.2*ampy])

bnm = interp(y,InterpFac,10,cutoff);

subplot(212)
stem(bnm,'bo');
ampbnm = max(abs(bnm));
ylabel(['Amplitude'])
xlabel(['(b)  Upsampled Seq. Using Interp'])
axis([0,InterpFac*SR,-1.2*ampbnm,1.2*ampbnm])

figure(35)

x = 0:1:SR-1;
xx = 0:1:SR*InterpFac-1;

subplot(311)
abfty = abs(fft(y,SR));
stem(xx(1,1:fix(length(abfty)/2)),abfty(1,1:fix(length(abfty)/2)),'bo');
xlabel(['(a)  Frequency'])
ylabel(['Mag'])
axis([0,fix(length(abfty)/2),0,1.2*max(abs(abfty)) ])

subplot(312)
abftnewy = abs(fft(Newy,SR*InterpFac));
stem(xx(1,1:fix(length(abftnewy)/2)),abftnewy(1,1:fix(length(abftnewy)/2)),'bo');
xlabel(['(b)  Frequency'])
ylabel(['Mag'])
axis([0 fix(length(abftnewy)/2) 0 1.2*max(abs(abftnewy))])

subplot(313)
abftOut = abs(fft(OutputInterp,SR*InterpFac));
stem(xx(1,1:fix(length(abftOut)/2)),abftOut(1,1:fix(length(abftOut)/2)),'bo');
xlabel(['(c)  Frequency'])
ylabel(['Mag'])
axis([0,fix(length(abftOut)/2),0,1.2*max(abs(abftOut))])
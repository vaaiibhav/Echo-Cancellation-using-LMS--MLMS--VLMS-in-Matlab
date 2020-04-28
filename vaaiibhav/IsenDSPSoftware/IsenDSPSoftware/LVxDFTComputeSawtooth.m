function LVxDFTComputeSawtooth
% This program demonstrates the computation of the Discrete Fourier
% Transform (DFT) of a Sawtooth wave using a direct implementation.
% as the input with unknown frequency content. 
%
% Press any key to compute the next bin.
% LVxDFTComputeSawtooth
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
SR = 32;
x = 0:1:SR-1;
t = 0:1/SR:1-1/SR;
DFTrealbin(1,1:SR) = 0;
DFTimagbin(1,1:SR) = 0;
Waveform(1,1:SR) = 0;

for n = 1:1:SR/2
   Waveform = Waveform + (1/n)*sin(2*pi*n*t);
end

A = 1/SR;
testvec = Waveform;
maxplot = 20/32;

figure(28);
clf

for harmno = 1:1:SR
testimag = -sin(2*pi*t*(harmno-1));
testreal = cos(2*pi*t*(harmno-1));
DFTrealbin(harmno) = A*sum(testvec.*testreal);
DFTimagbin(harmno) = A*sum(testvec.*testimag);

subplot(321)
stem(t*32,testreal,'bo');
ylabel(['Amp'])
xlabel(['(a) Sample'])
axis([-1 32 -1.5 1.5])

subplot(323)
stem(t*32,testimag,'bo');
ylabel(['Amp'])
xlabel(['(c) Sample'])
axis([-1 32 -1.5 1.5])

subplot(325)
stem(t*32,testvec,'bo');
ylabel(['Amp'])
xlabel(['(e) Sample'])
axis([-1 32  -1.2*abs(min(testvec)) 1.2*max(testvec)])

subplot(322)
stem(x,DFTrealbin,'bo');
ylabel(['Real'])
xlabel(['(b) Bin'])
axis([-1 32 -maxplot maxplot])

subplot(324)
stem(x,DFTimagbin,'bo');
ylabel(['Imag'])
xlabel(['(d) Bin'])
axis([-1 32 -maxplot maxplot])

subplot(326)
njk = (DFTrealbin.^2 + DFTimagbin.^2).^0.5;
stem(x,njk,'bo');
ylabel(['Mag'])
xlabel(['(f) Bin'])
axis([-1 32 0 maxplot])

if harmno==32
   return
else
   pause%(0.5) 
end

end
function LVxDFTComputeImpUnBal
% function DemoDFTComputeImpUnBal
% This program demonstrates the computation of the Discrete Fourier
% Transform (DFT) of an imipulse fcn over N samples using a direct
% implementation.
% The call:
% LVxDFTComputeImpUnBal
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
t = 0:1/32:1-1/32;

testvec = [1,  zeros(1,31)];

A = 1/32;
DFTrealbin = zeros(1,32);
DFTimagbin = zeros(1,32);
x = 1:1:32;

figure(27);

for harmno = 1:1:32
testimag = -sin(2*pi*t*(harmno-1));
testreal = cos(2*pi*t*(harmno-1));
DFTrealbin(1,harmno) = A*sum(testvec.*testreal);
DFTimagbin(1,harmno) = A*sum(testvec.*testimag);

subplot(321)
stem(t*32,testreal,'bo');
ylabel(['Amplitude'])
xlabel(['(a)  Test Cosine-Bin ',num2str(harmno-1),'; X-Axis = Samp No.'])
axis([-1 32 -1.3 1.3])

subplot(323)
stem(t*32,testimag,'bo');
ylabel(['Amplitude'])
xlabel(['(c)  Test Sine-Bin ',num2str(harmno-1),'; X-Axis = Samp No.'])
axis([-1 32 -1.3 1.3])

subplot(325)
stem(t*32,testvec,'bo');
ylabel(['Amplitude'])
xlabel(['(e) Test Signal; X-Axis = Samp No.'])
axis([-1 32 -1.3 1.3])

subplot(322)
stem(x-1,DFTrealbin,'bo');
ylabel(['Real(DFT)'])
xlabel(['(b)  Bin'])
axis([-1 32 -0.05 0.05])

subplot(324)
stem(x-1,DFTimagbin,'bo');
ylabel(['Imag(DFT)'])
xlabel(['(d)  Bin'])
axis([-1 32 -0.05 0.05])

subplot(326)
njk = (DFTrealbin.^2 + DFTimagbin.^2).^0.5;
stem((x-1),njk,'bo');
ylabel(['DFT Mag'])
xlabel(['(f)  Bin'])
axis([-1 32 0 0.05])

if harmno==32
   return
else
   pause
end

end
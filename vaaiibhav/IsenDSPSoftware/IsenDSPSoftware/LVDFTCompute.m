function LVDFTCompute(TestSignal)
% function LVDFTCompute(TestSignal)
%
% This program demonstrates the computation of the Discrete Fourier
% Transform (DFT)using a direct implementation.
%
% TestSignal = 1 uses a series of ones followed by a
% series of negative ones to form a finite length square wave;
% TestSignal = 0 synthesizes a partial square wave by summing
% odd harmonics only up to half the sequence length of 32.
% The DFT is computed bin by bin by pressing any key to compute the next bin.
%
% The two possible calls are: 
%
% LVDFTCompute(0)
% LVDFTCompute(1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

t = 0:1/32:1-1/32;
if TestSignal==0
Signal = sin(2*pi*t)+ 0.3333*sin(2*pi*t*3)+ 0.2*sin(2*pi*t*5)+...
 0.1428*sin(2*pi*t*7)+ 0.1111*sin(2*pi*t*9)+ 0.0909*sin(2*pi*t*11)+...
 0.0769*sin(2*pi*t*13) + 0.0666*sin(2*pi*t*15);
elseif TestSignal==1
   Signal = [ones(1,16),  -ones(1,16)];
else
   Signal = [ones(1,16),  -ones(1,16)];   
end

limharm = 32;
A = 1/32;
F = zeros(1,limharm);
x = 1:1:limharm;

figure(21);
for k = 1:1:limharm
TestCExp = exp(-j*2*pi*t*(k-1));
testimag = -sin(2*pi*t*(k-1));  % for display
testreal = cos(2*pi*t*(k-1));   % for display

F(k) = A*sum(Signal.*TestCExp);

subplot(321)
stem(t*32,testreal,'bo');
ylabel(['Amplitude'])
xlabel(['(a)  Test Cosine-Bin ',num2str(k-1)])
axis([0,length(t),-1.3,1.3])

subplot(323)
stem(t*32,testimag,'bo');
ylabel(['Amplitude'])
xlabel(['(c)  Test Sine-Bin ',num2str(k-1)])
axis([0, length(t),-1.3, 1.3])

subplot(325)
stem(t*32,Signal,'bo');
xlabel(['(e) Signal'])
ylabel(['Amplitude'])
axis([0, length(t), -1.3, 1.3])

subplot(322)
stem(x-1,real(F),'bo');
ylabel(['Real(F[k])'])
xlabel(['(b)  Bin (k)'])
axis([0, limharm, -1.2, 1.2])

subplot(324)
stem(x-1,imag(F),'bo');
ylabel(['Imag(F[k])'])
xlabel(['(d)  Bin (k)'])
axis([0, limharm, -1.2, 1.2])

subplot(326)
njk = abs(F);
stem((x-1),njk,'bo');
axis([0 limharm 0 1.2])
xlabel(['(f)  Bin (k)'])
ylabel(['Mag(F[k])'])
%hold on
if k==limharm
   return
else
   pause
end

end




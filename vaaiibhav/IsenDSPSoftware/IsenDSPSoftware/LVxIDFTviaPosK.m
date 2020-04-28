function LVxIDFTviaPosK(TestSig)
% function LVxIDFTviaPosK(TestSig)
% Test signal must be real
% Computes the dft of TestSig, then computes the IDFT
% two ways, using a formula (see text) for the IDFT (for real signal sequences only,
% assumes DFT bins are complex conjugates of each other except for Bins 0
% and N/2;
% TestSig is also reconstructed via standard IDFT using the
% function ifft, and a figure with three plots is generated, the first plot
% showing the TestSig, the second one showing the reconstruction via
% the text formula for IDFT for real signals, and the third is TestSig
% reconstructed via the function ifft.
% Test calls:
% LVxIDFTviaPosK(ones(1,8))
% LVxIDFTviaPosK(ones(1,9))
% LVxIDFTviaPosK(randn(1,9))
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
N = length(TestSig);

ft = fft(TestSig);

if rem(N,2)==0  % N is even
    k = 0:1:N/2;
else
    k = 0:1:(N-1)/2;
end
n = 0:1:N-1;
Nov2 = N/2;
x = 0;
for hCtr = k
    if hCtr == 0
  %      x = x + ft(hCtr+1)*cos(2*pi*n*0/N);
        x = x + ft(hCtr+1);
    elseif hCtr==Nov2
 %       x = x + ft(Nov2+1)*cos(2*pi*n*(Nov2)/N);
 %       x = x + ft(Nov2+1)*cos(n*pi);
        x = x + ft(Nov2+1)*(-1).^n;
    else
        x = x + 2*( real(ft(hCtr+1))*cos(2*pi*n*hCtr/N) - imag(ft(hCtr+1))*sin(2*pi*n*hCtr/N));
    end
end
x = (1/N)*x;

figure(777)

subplot(311)
stem(n,TestSig)
xlabel('(a) n (TestSig)')
ylabel('Amplitude')
axis([0, length(TestSig),min(TestSig)-1,max(TestSig)+1])

subplot(312)
stem(n,x)
xlabel('(b) n (Reconstruction of TestSig from DFT Via Text Formula)')
ylabel('Amplitude')
axis([0, length(TestSig),min(TestSig)-1,max(TestSig)+1])

subplot(313)
xx = real(ifft(ft));
stem(n,xx)
xlabel('(c) n (Reconstruction of TestSig from DFT Via function ifft)')
ylabel('Amplitude')
axis([0, length(TestSig),min(TestSig)-1,max(TestSig)+1])
    
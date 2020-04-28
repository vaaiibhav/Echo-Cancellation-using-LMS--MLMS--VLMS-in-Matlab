function [seqCos,seqSin] = LVGenFreq(M,k,N)
% [seqCos,seqSin] = LVGenFreq(1,7.5,73)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
n = [0:1:N-1]; arg = 2*pi*k/N;
mags = (M.^n); maxmags = max(abs(mags));
W2n = mags.*exp(j*arg).^n;
seqCos = real(W2n); seqSin = imag(W2n); 
figure(65); subplot(211); stem(n,seqCos);
xlabel('(a) Sample'); ylabel('Amplitude');
axis([0,N-1,-1.1*maxmags,1.1*maxmags])
subplot(212); stem(n,seqSin);
xlabel('(b) Sample'); ylabel('Amplitude')
axis([0,N-1,-1.1*maxmags,1.1*maxmags])
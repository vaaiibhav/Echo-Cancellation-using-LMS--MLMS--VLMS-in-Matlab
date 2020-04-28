function LVxOrthogAudio(DecodePhi)
% Creates sine and cosine carriers of frequency greater than 16 kHz, modulates each by a corresponding 
% audio file drwatsonSR8K.wav or whoknowsSR8K.wav, takes the difference to create a signal 
% transmission signal, then decodes the transmission signal to produce one or the other of the encoded 
% audio signals, according to the value of the variable DecodePhi in
% the function call.
% DecodePhi is the phase angle of the decoding carrier, and should be
% between 0 and pi/2; 0 will decode the sine carrier's audio, pi/2 will
% decode the cosine carrier's audio, and numbers in between 0 and pi/2 will
% cause a proportional mixture of the two audio audio signals to be
% decoded.
% LVxOrthogAudio(0)
% LVxOrthogAudio(pi/2)
% LVxOrthogAudio(pi/4)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

[A,Fs,NBITS] = wavread('drwatsonSR8K.wav');
[B,Fs,NBITS] = wavread('whoknowsSR8K.wav');
A = A'; B = B';
themin = min([length(A), length(B)]);
A = A(1,1:themin); B = B(1,1:themin);    
N = 10; F=1;   
n = 0:1:N-1;
C1 = cos(2*pi*n*F/N);
C2 = sin(2*pi*n*F/N);
C1Mat = (C1')*(ones(1,themin));
C2Mat = (C2')*(ones(1,themin));
AMat = ones(N,1)*A; BMat = ones(N,1)*B;     
S1 = C1Mat.*(AMat); S1 = S1(:);
S2 = C2Mat.*(BMat); S2 = S2(:);
S = S1 - S2;
% must break S into one cycle frames
SigMat = reshape(S,N,fix(length(S)/N));
DecodeMat = (cos(2*pi*n*F/N + DecodePhi)')*(ones(1,themin));
Sr = (2/N)*sum(SigMat.*DecodeMat);

figure(333)

subplot(321)
plot(A)
xlabel('(a) n')

subplot(322)
plot(S1)
xlabel('(b) n')

subplot(323)
plot(B)
xlabel('(c) n')

subplot(324)
plot(S2)
xlabel('(d) n')

subplot(325)
plot(S)
xlabel('(e) n')

subplot(326)
plot(Sr)
xlabel('(f) n')

Sr = 0.99*Sr/max(abs(Sr));
sound(Sr,8000)

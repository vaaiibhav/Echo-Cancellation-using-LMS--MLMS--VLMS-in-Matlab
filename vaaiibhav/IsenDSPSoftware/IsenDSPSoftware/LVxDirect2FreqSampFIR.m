function [CFsB,CFsA,BFs,AFs] = LVxDirect2FreqSampFIR(Imp)
% function [CFsB,CFsA,BFs,AFs] = LVxDirect2FreqSampFIR(Imp)
% Receives an FIR impulse response Imp and generates the Frequency Sampling
% Coefficients BFs,and AFs with the comb filter section coefficients as CFsB 
% and CFsA.
% Sample calls:
% [CFsB,CFsA,BFs,AFs] = LVxDirect2FreqSampFIR([1,1,1,1])
% [CFsB,CFsA,BFs,AFs] = LVxDirect2FreqSampFIR([1,-1,1,-1])
% [CFsB,CFsA,BFs,AFs] = LVxDirect2FreqSampFIR([1,0,0,1])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
N = length(Imp)
Hz = fft(Imp)
angHk = angle(Hz);

H0 = Hz(1)
if rem(N,2)==0 % even length, will have N/2 sample
    HN2 = Hz(N/2+1)
    M = N/2 - 1;
    BFs = zeros(M+2,2);
    AFs = zeros(M+2,3);
    BFs(M+1,:) = [H0, 0];
    AFs(M+1,:) = [1,-1,0];
    BFs(M+2,:) = [HN2,0];
    AFs(M+2,:) = [1,1,0];
else
    M = (N-1)/2;
    BFs = zeros(M+1,2)
    AFs = zeros(M+1,3);
    BFs(M+1,:) = [H0,0];  
    AFs(M+1,:) = [1,-1,0];
end

for k = 1:1:M;
    BFs(k,:) = 2*abs(Hz(k+1))*[ cos(angHk(k+1)),-cos( angHk(k+1)-2*pi*k/N ) ];
    AFs(k,:) = [1,-2*cos(2*pi*k/N),1];
end

CFsB = [1,zeros(1,N-1),-1];
CFsA = [N];

    
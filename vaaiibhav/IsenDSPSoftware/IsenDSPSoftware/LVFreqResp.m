function LVFreqResp(tstSig, NoFreqs)
% LVFreqResp([1.9, 0, -0.9, 0, 1.9, 0, -3.2], 500)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
NyqLim = length(tstSig)/2; LTS = length(tstSig);
t = [0:1:(LTS-1)]/(LTS); FR = [];
frVec = 0:NyqLim/(NoFreqs-1):NyqLim;
    for Freq = frVec     % FR via loop
    testCorr = cos(2*pi*t*Freq) - j*sin(2*pi*t*Freq);
    FR = [FR, sum(tstSig.*testCorr)]; 
    end
% FR via vectorized operation
% FR = exp(-j*(((2*pi*t)'*frVec)'))*(tstSig'); 
figure(9)
xvec = frVec/(frVec(length(frVec)));
plot(xvec,abs(FR));
xlabel('Normalized Frequency')
ylabel('Magnitude')
axis([0,inf,0,inf])
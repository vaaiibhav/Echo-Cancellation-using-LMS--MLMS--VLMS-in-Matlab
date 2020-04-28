function [FR] = LVxFreqRespND(tstSig, LenCorr)
% FR = LVxFreqRespND([ones(1,32)], 128)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
LTS = length(tstSig);
if LenCorr>LTS
    tstSig = [tstSig,zeros(1,LenCorr-LTS)];
elseif LenCorr < LTS
    error('LenCorr must be at least equal to the length of tstSig')
end
ultLen = length(tstSig);
t = [0:1:(ultLen-1)]/ultLen; 
if rem(LenCorr,2)==0   
    Lim = LenCorr/2;
else
    Lim = (LenCorr-1)/2;
end
frVec = 0:1:Lim;
FR = exp(-j*(((2*pi*t)'*frVec)'))*(tstSig'); 


function [OutputSeq] = LVxDFTorIDFT(InputSeq, DftOrIdftFLAG)
% A is a row vector, which is either a time domain sequence to have the DFT performed on it or
% a DFT to be turned into a time domain sequence via the IDFT, which is performed with 
% the same DFT algorithm
% Pass DftOrIdftFLAG as 0 to perform the DFT or 1 to perform the IDFT
% Sample Calls:
%  [OutputSeq] = LVxDFTorIDFT([ones(1,8)], 0)
%  [OutputSeq] = LVxDFTorIDFT([ones(1,8)], 1)
%  [OutputSeq] = LVxDFTorIDFT([8, zeros(1,7)], 0)
%  [OutputSeq] = LVxDFTorIDFT([8, zeros(1,7)], 1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

N = length(InputSeq);
n = 0:1:N-1;
k = n;

OutputSeq = zeros(1,N);

if DftOrIdftFLAG==0 % dft
else
    InputSeq = conj(InputSeq);
end

for k = 1:1:N;
    OutputSeq(1,k:k) = sum(InputSeq.*exp(-j*2*pi*n*(k-1)/N));
end

if DftOrIdftFLAG==0 % dft
else
    OutputSeq = (1/N)*conj(OutputSeq);
end



function [DFTGBin,MagDFTGBin] = LVxDFTViaGoertzelBin(Bin,Signal)
% [DFTGBin,MagDFTGBin] = LVxDFTViaGoertzelBin(13,randn(1,31))
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
x = Signal;
N=length(x); 
GB = filter(1,[1  -2*(cos(2*pi*Bin/N))  1],[x 0]);
GrtzlBin = GB(N+1) - exp(-j*2*pi*Bin/N)*GB(N);
MagSqGBin = GB(N+1)^2 - 2*cos(2*pi*Bin/N)*GB(N+1)*GB(N) + GB(N)^2;

MagDFTGBin = MagSqGBin^0.5;
DFTGBin = GrtzlBin;

%ft = fft(x); 
%FFTBin = ft(Bin+1);
%MagFFTBin = abs(FFTBin);


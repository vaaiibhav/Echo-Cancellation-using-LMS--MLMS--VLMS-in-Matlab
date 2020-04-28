function LVTypeIVHilbertViaIDFT(ZeroBin,PosBins,NOver2Bin,NegBins)
%
% LVTypeIVHilbertViaIDFT(0,[-j*ones(1,9)],[-j],[j*ones(1,9)])
% LVTypeIVHilbertViaIDFT(0,[-j*ones(1,29)],[-j],[j*ones(1,29)])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Bins = [ZeroBin,PosBins,NOver2Bin,NegBins];
L = length(Bins), if ~(rem(L,2)==0);
error('Total Bin Count is not an even number'); end
M = (L-1)/2; kZeroPosBns = 0:1:ceil(M); 
kNegBns = ceil(M)+1:1:L-1;
angVec = [-2*pi*kZeroPosBns*M/L,2*pi*(L-kNegBns)*M/L];
PhaseFac = exp(j*angVec); NetBins = Bins.*PhaseFac; 
Imp = real(ifft(NetBins)); figure(88); stem(Imp)
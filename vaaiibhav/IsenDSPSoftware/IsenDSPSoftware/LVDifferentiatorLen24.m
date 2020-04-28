function LVDifferentiatorLen24
% LVDifferentiatorLen24
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Bins = pi*( [ 0, j*(1:1:11)*2/24, j*12*2/24, j*(-11:1:-1)*2/24 ] )
L = length(Bins); M = (L-1)/2;
kZeroPosBns = 0:1:ceil(M); kNegBns = ceil(M)+1:1:L-1;
angVec = [-2*pi*kZeroPosBns*M/L,2*pi*(L-kNegBns)*M/L];
PhaseFac = exp(j*angVec); NetBins = Bins.*PhaseFac; 
Imp = (1/pi)*real(ifft(NetBins))
figure(98); stem(Imp); xlabel('Sample'); ylabel('Amplitude')
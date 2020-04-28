function LVTypeIVHilbertViaIDFTSimple(L)
% LVTypeIVHilbertViaIDFTSimple(60)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(rem(L,2)==0)
L = L + 1; end
MM = L/2 - 1;
Bins = [ [j*ones(1,MM)],[0],[-j*ones(1,MM)],[-j] ];
M = (L-1)/2; k = -L/2+1:1:L/2; LNB = (L-2)/2; 
PhaseFac = exp(-j*2*pi*k*M/L);
NetBns = Bins.*PhaseFac; 
NetBns = [NetBns(1,LNB+1:L),NetBns(1,1:LNB)];
Imp = real(ifft(NetBns)); figure(89); stem(Imp)
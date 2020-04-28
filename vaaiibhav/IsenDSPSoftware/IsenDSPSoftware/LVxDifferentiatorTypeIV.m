function ImpResp = LVxDifferentiatorTypeIV(L) 
% LVxDifferentiatorTypeIV(24)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(rem(L,2)==0)
    L = L + 1;
end
Bins = pi*( [ 0, j*(1:1:L/2-1)*2/L, j*(L/2)*2/L, j*(-(L/2-1):1:-1)*2/L ] );
L = length(Bins); M = (L-1)/2;
kZeroPosBns = 0:1:ceil(M); kNegBns = ceil(M)+1:1:L-1;
angVec = [-2*pi*kZeroPosBns*M/L,2*pi*(L-kNegBns)*M/L];
PhaseFac = exp(j*angVec); 
NetBins = Bins.*PhaseFac; ImpResp = (1/pi)*real(ifft(NetBins));
fr = fft(ImpResp,2048); fr = abs(fr(1,1:1024));
figure(98); subplot(211); stem(ImpResp);
subplot(212); plot(fr)
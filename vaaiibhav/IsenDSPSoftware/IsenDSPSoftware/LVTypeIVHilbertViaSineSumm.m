function LVTypeIVHilbertViaSineSumm(L)
% LVTypeIVHilbertViaSineSumm(60)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(rem(L,2)==0); L = L + 1; end;
Ak0 = 0; AkPos = [0.8,ones(1,L/2-2)]; AkLOver2 = 1;
M = (L-1)/2; Ak = [Ak0,AkPos,AkLOver2];
limK = L/2; WF = zeros(1,L); n = 0:1:L-1;
for k = 1:1:limK;
if (k==limK); C = 1; else; C = 2; end; 
WF = WF + C*Ak(k+1)*sin(2*pi*(n-M)*k/L);
end; WF = WF/L; figure(99); stem(WF)
function LVCheby1MagSquared(Ep,N)  
% LVCheby1MagSquared(0.4,5)  
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
inc = 0.01; xLo = 0:inc:1;
xHi = 1+inc:inc:3; TnLo = cos(N*acos(xLo));
TnHi = cosh(N*acosh(xHi)); T = [TnLo TnHi];
MagHSq = 1./(1 + Ep^2*(T.^2));
figure; xplot = [xLo, xHi]; plot(xplot,MagHSq); 
xlabel('Norm Freq (\Omega/\Omegac)'); ylabel('Mag Squared')
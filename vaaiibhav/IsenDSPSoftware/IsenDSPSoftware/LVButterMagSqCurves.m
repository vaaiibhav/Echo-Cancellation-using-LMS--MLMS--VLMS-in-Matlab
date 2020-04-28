function LVButterMagSqCurves
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
figure(79); 
NormFrq = 0:0.01:2;
for N = 1:1:12; hold on
H = 1./(1 + (NormFrq).^(2*N)); 
plot(NormFrq,H,'k--'); end; hold off
xlabel('Frequency, Units of \Omega_C')
ylabel('Magnitude')


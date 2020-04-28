function LVxDigEllipBPFViaMS(Rp,As,ws1,wp1,wp2,ws2)
% Receives the values for maximum passband ripple Rp, minimum stopband
% attenuation As, and bandpass filter band edges ws1,wp1,wp2, and ws2, and
% designs a digital bandpass filter using the script LVxDesignDigEllipBPF
% and the built-in function ellip. Plots the magnitude and phase responses
% of both designed filters; the realized values for Rp and As are also
% computed for both filters.
% Test calls:

% LVxDigEllipBPFViaMS(1,75,0.4,0.475,0.65,0.775)
% LVxDigEllipBPFViaMS(1,45,0.4,0.475,0.65,0.775)
% LVxDigEllipBPFViaMS(1,60,0.45,0.55,0.7,0.83)
% LVxDigEllipBPFViaMS(1,75,0.4,0.475,0.65,0.775)
% LVxDigEllipBPFViaMS(1,75,0.1,0.12,0.4,0.5)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
[b,a,G,NetRp,NetAs] = LVxDesignDigEllipBPF(Rp,As,ws1,wp1,wp2,ws2)

%NN = length(roots(a))/2

NN = LVxOrdEllipBPF(Rp,As,ws1,wp1,wp2,ws2)

[bMS,aMS] = ellip(NN,Rp,As,[wp1,wp2]);

NumMSEllipBPFPoles = length(roots(aMS))

Hz = LVzFr(bMS,aMS,2000,46);
[NetRpMS,NetAs1MS,NetAs2MS] = LVxRealizedFiltParamBPF(Hz,ws1,wp1,wp2,ws2,1)
function [zB,zA,G,NetRp,NetAs] = LVxDesignDigitalButterLPF(Rp,As,wp,ws)
% function [zB,zA,G,NetRp,NetAs] = LVxDesignDigitalButterLPF(Rp,As,wp,ws)
% Designs a digital Butterworth LPF having no more than Rp dB of ripple at
% passband edge wp, and at least As dB attenuation at stopband edge ws.
% The filter coefficients are returned as numerator coefficients zB,
% denominator coefficeints zA, and gain G. The script also plots the
% magnitude and phase responses of both the digital filter, and the analog
% prototype used to compute the digital filter coefficients using the
% Bilinear transform.
% Sample calls:
% [zB,zA,G,NetRp,NetAs] = LVxDesignDigitalButterLPF(1,40,0.5*pi,0.6*pi)
% [zB,zA,G,NetRp,NetAs] = LVxDesignDigitalButterLPF(0.3,60,0.25*pi,0.325*pi)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

Fs = 1; T = 1/Fs; OmP = (2/T)*tan(wp/2);
OmS = (2/T)*tan(ws/2); E = sqrt(10^(Rp/10)-1);
A = 10^(As/20); 

[Z,P,K] = LVDesignButterworth(OmP,OmS,Rp,As)

OmC = OmP; % OmT = OmS/OmP; 
%g = sqrt((A^2-1)/(E^2)); 
%NN = ceil(log10(g + sqrt(g^2 - 1))/log10(OmT + sqrt(OmT^2-1)))
%v0 = asinh(1/E)/NN; k = -(NN-1):2:(NN-1);
%P = OmC*(-sinh(v0)*cos(k*pi/(2*NN)) + ...
%    j*cosh(v0)*sin(k*pi/(2*NN)));
%NetK = prod(abs(P));
%if rem(NN,2)==0 
%NetK = NetK/sqrt(1 + E^2); end


sB = K; sA = poly(P); Num = K; Den = 1; 
for Ctr = 1:1:length(P); D = [1 1]; 
N = [(2*Fs-P(Ctr)) -(2*Fs+P(Ctr))];
Num = conv(Num,N); Den = conv(Den,D); 
end; 
zB = real(Den); zA = real(Num);
b0 = 1/zB(1); a0 = 1/zA(1);
zA = zA*a0; zB = zB*b0;
G = abs(polyval(zA,exp(j*0))./polyval(zB,exp(j*0)));
 
zB = G*zB; 
LVsFRzFr(sB,sA,OmC,zB,zA,44); 
Hz = LVzFr(zB,zA,1000,45);
[NetRp,NetAs] = LVzRealizedFiltParams(Hz,wp,ws,pi);


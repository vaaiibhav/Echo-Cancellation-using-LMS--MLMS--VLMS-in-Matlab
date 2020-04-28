function LVsCheby1LPF2zCheby1LPF(Rp,As,wp,ws)
% LVsCheby1LPF2zCheby1LPF(1,40,0.5*pi,0.6*pi)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Fs = 1; T = 1/Fs; OmP = (2/T)*tan(wp/2);
OmS = (2/T)*tan(ws/2); E = sqrt(10^(Rp/10)-1);
A = 10^(As/20); OmC = OmP; OmT = OmS/OmP; 
g = sqrt((A^2-1)/(E^2)); 
NN = ceil(log10(g + sqrt(g^2 - 1))/log10(OmT + sqrt(OmT^2-1)))
v0 = asinh(1/E)/NN; k = -(NN-1):2:(NN-1);
P = OmC*(-sinh(v0)*cos(k*pi/(2*NN)) + ...
    j*cosh(v0)*sin(k*pi/(2*NN)));
NetK = prod(abs(P));
if rem(NN,2)==0 
NetK = NetK/sqrt(1 + E^2); end
sB = NetK; sA = poly(P); Num = NetK; Den = 1; 
for Ctr = 1:1:length(P); D = [1 1]; 
N = [(2*Fs-P(Ctr)) -(2*Fs+P(Ctr))];
Num = conv(Num,N); Den = conv(Den,D); 
end; 
zB = real(Den); zA = real(Num);
b0 = 1/zB(1); a0 = 1/zA(1);
zA = zA*a0; zB = zB*b0;
G = abs(polyval(zA,exp(j*0))./polyval(zB,exp(j*0)));
if rem(NN,2)==0; G = G/sqrt(1 + E^2); end 
zB = G*zB; 
LVsFRzFr(sB,sA,OmC,zB,zA,44); Hz = LVzFr(zB,zA,1000,45);
[NetRp,NetAs] = LVzRealizedFiltParamLPF(Hz,wp,ws,pi)


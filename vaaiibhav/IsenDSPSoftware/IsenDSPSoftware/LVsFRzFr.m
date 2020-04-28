function LVsFRzFr(sB,sA,OmegaC,zB,zA,FigNo)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
LnLm = 8*OmegaC;
rootrat = ((10^6)/(LnLm))^(0.001);
Sargs = 0:LnLm/1000:LnLm-LnLm/1000;
yy = (LnLm)*(rootrat.^[1:1:1000]);
Sargs = [Sargs, yy]; s = j*Sargs;
Hs = polyval(sB,s)./polyval(sA,s);
Zargs = 0:0.01:pi; z = exp(j*Zargs);
Hz = polyval(zB,z)./polyval(zA,z);
figure(FigNo); subplot(211); semilogx(Sargs,abs(Hs))
xlabel('(a) Freq, Radians/s'); ylabel('Magnitude')
%axis([0,max(Sargs),0,1.1])
subplot(212); plot(Zargs/pi,abs(Hz)); grid on;
xlabel('(b) Freq, Units of \pi'); ylabel('Magnitude')
axis([0,1,0,1.1])
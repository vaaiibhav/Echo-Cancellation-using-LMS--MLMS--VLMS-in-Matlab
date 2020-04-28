function LVsFRzFrLog(sB,sA,OmegaC,zB,zA,FigNo)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
LnLm = 8*OmegaC;
Sargs = 0:0.01:LnLm;
Sargs = [Sargs]; s = j*Sargs;
Hs = polyval(sB,s)./polyval(sA,s);
Zargs = 0:0.01:pi; z = exp(j*Zargs);
Hz = polyval(zB,z)./polyval(zA,z);

figure(FigNo); subplot(211); plot(Sargs,20*log10(abs(Hs+eps)))
xlabel('(a) Freq, Radians/s'); ylabel('Magnitude, dB');
axis([0 LnLm -100 5]); 
subplot(212); 
ploty = 20*log10(abs(Hz)+eps);
plot(Zargs/pi,ploty); grid on;
xlabel('(b) Freq, Units of \pi'); ylabel('Magnitude, dB');
axis([10^(-2) inf -100 max(ploty)+10])
function LVSsbDsb(Fc,Fa,Phi1,Phi2,Sb,SR)
% LVSsbDsb(100,20,pi/2,pi/6,0,1000)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
phi1 = pi/2; phi2 = pi/6; t = 0:1/SR:1-1/SR; 
argsC = 2*pi*t + phi1; argsS = 2*pi*t + phi2; 
dsb = cos(argsC*Fc).*cos(argsS*Fa); 
figure(95); subplot(2,1,1); 
plot(2*t,abs(fft(dsb))); xlabel('Freq, Units of \pi')
ylabel('Mag, DSB Signal')
if Sb==0; ssb = dsb + sin(argsC*Fc).*sin(argsS*Fa);
else; ssb = dsb - sin(argsC*Fc).*sin(argsS*Fa); end
subplot(2,1,2); plot(2*t,abs(fft(ssb)))
xlabel('Freq, Units of \pi'); ylabel('Mag, SSB Signal')
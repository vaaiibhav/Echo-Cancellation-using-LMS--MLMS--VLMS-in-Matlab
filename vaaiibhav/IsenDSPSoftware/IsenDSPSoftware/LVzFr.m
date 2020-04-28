function Hz = LVzFr(zB,zA,NumPts,FigNo)
% Computes and displays on Figure FigNo the magnitude and phase response of a digital
% filter having numerator coefficients zB and denomiantor coefficients zA,
% for NumPts frequency points. The complex frequency response vector Hz is
% returned by the function and can be supplied to another script if desired
% to compute realized filter parameters.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Zargs = 0:pi/NumPts:pi; z = exp(j*Zargs);
Hz = polyval(zB,z)./polyval(zA,z);
figure(FigNo); subplot(211)
plot(Zargs/pi,20*log10(abs(Hz)+eps)); grid
xlabel('(a) Freq, Units of \pi Radians'); 
ylabel('Magnitude, dB'); axis([0,1,-100,10])
subplot(212)
plot(Zargs/pi,unwrap(angle(Hz))); grid
xlabel('(b) Freq, Units of \pi Radians'); ylabel('Radians')
axis([0,1,-inf,inf])
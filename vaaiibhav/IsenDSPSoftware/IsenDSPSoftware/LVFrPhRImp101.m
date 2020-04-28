function LVFrPhRImp101
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
incF = 2*pi/1024; argF = -pi+incF:incF:pi;
Hr = 2*cos(argF); Hph = -argF;
figure(8); subplot(221); plot(argF/pi, Hr)
subplot(222); plot(argF/pi,Hph)
H = fft([1 0 1],1024); Hdtft = fftshift(H); subplot(223); 
xvec = [-511:1:512]/512; plot(xvec, abs(Hdtft)); 
subplot(224); plot(xvec, unwrap(angle(Hdtft)))
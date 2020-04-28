function LVxDesignDigCheby1BPFViaImpInv(Ws1,Wp1,Wp2,Ws2,Rp,As,Fs)
% Designs a digital bandpass Chebyshev Type I filter using
% impulse invariance. The four bandpass filter band limits are
% Ws1,Wp1,Wp2,and Ws2. The maximum allowable passband ripple in positive dB 
% is Rp and the minimum stopband attenuation in positive dB is As.
% Test calls:
% LVxDesignDigCheby1BPFViaImpInv(0.2*pi,0.3*pi,0.6*pi,0.8*pi,0.5,40,0.05)
% LVxDesignDigCheby1BPFViaImpInv(0.3*pi,0.38*pi,0.6*pi,0.72*pi,0.5,50,0.05)
% LVxDesignDigCheby1BPFViaImpInv(0.2*pi,0.35*pi,0.6*pi,0.9*pi,1,60,0.1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
OmS1=Ws1*Fs;
OmP1=Wp1*Fs;
OmP2=Wp2*Fs;
OmS2=Ws2*Fs;

[Z,P,K] = LVDesignCheby1BPF(OmS1,OmP1,OmP2,OmS2,Rp,As);   
b = real(K*poly(Z)), a = real(poly(P)),
[BZ,AZ] = impinvar(b,a,Fs);
LVsFRzFrLog(b,a,OmP2,BZ,AZ,97)
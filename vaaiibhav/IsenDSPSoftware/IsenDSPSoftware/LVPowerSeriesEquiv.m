function LVPowerSeriesEquiv(M,k,N)
% LVPowerSeriesEquiv(0.9,3,64)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
n = [0:1:N-1]; arg = 2*pi*k/N;
mags = (M.^n); maxmags = max(mags);
W2n = mags.*exp(j*arg).^n;
rightS = mags.*(cos(n*arg) + j*sin(n*arg));
figure(19); clf; hold on;
plot(real(W2n),imag(W2n),'bo');
plot(real(rightS),imag(rightS),'rx');
grid on; xlabel('Real'); ylabel('Imag')
axis([-maxmags,maxmags,-maxmags,maxmags])

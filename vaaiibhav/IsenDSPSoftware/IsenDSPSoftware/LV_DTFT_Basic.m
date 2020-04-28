function LV_DTFT_Basic(x,M,R)
% LV_DTFT_Basic([1 0 1],300,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
N = length(x); W = exp(-j*R*pi/M); k = 0:1:M-1;
n = 0:1:N-1; dMat = W.^(n'*k); d = x*dMat; figure(9)
subplot(2,2,1); plot(R*[0:1:M-1]/M,abs(d)); 
grid; xlabel('Norm Freq'); ylabel('Mag')
subplot(2,2,2); plot(R*[0:1:M-1]/M,angle(d))
grid; xlabel('Norm Freq'); ylabel('Radians')
subplot(2,2,3); plot(R*[0:1:M-1]/M,real(d)); 
grid; xlabel('Norm Freq'); ylabel('Real')
subplot(2,2,4); plot(R*[0:1:M-1]/M,imag(d))
grid; xlabel('Norm Freq'); ylabel('Imag')
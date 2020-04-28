function LVButterPoles(N,OmegaC)
% LVButterPoles(3,2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
k = 0:1:2*N-1; 
P = OmegaC*exp(j*(pi/(2*N))*(2*k+N+1))
NetP = P(find(real(P)<0))
figure(90); clf; hold on; args = 0:0.02:2*pi;
plot(real(NetP),imag(NetP),'bx');
xlabel('Real')
cnums = OmegaC*exp(j*args);
plot(real(cnums),imag(cnums),':')
ylabel('Imaginary')
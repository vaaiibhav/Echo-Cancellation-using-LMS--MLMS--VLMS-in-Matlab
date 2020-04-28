function LVCheby1(N,OmC,Epsilon) 
% LVCheby1(5,1,0.5)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Alpha = 1/Epsilon + sqrt(1 + 1/(Epsilon^2));
B = 0.5*( Alpha^(1/N) + (1/Alpha)^(1/N) );
A = 0.5*( Alpha^(1/N) - (1/Alpha)^(1/N) );
k = 0:1:N-1; arg = pi/2 + pi*(2*k + 1)/(2*N);
SigK = A*cos(arg); OmK = B*sin(arg);
P = OmC*(SigK + j*OmK); K = prod(abs(P));
if rem(N,2)==0 % N even
K = 1/sqrt(1 + Epsilon^2)*K; end;
H = LVsFreqResp(K,poly(P),2*OmC,8);
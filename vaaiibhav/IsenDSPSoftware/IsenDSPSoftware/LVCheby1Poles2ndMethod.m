function p = LVCheby1Poles2ndMethod(N,OmC,Epsilon)
% p = LVCheby1Poles2ndMethod(5,2,0.5)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
k = -(N-1):2:N-1; r = k*pi/(2*N);
v0 = asinh(1/Epsilon)/N;
p = OmC*(-sinh(v0)*cos(r) + j*cosh(v0)*sin(r));
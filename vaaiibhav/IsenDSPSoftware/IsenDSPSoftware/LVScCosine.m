function LVScCosine(A,N,F)
% LVScCosine(1,128,3)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
t = [0:1:N-1]/N; x = [A*cos(2*pi*F*t)];y = 2*x; 
subplot(2,1,1); stem(x); subplot(2,1,2); stem(y)
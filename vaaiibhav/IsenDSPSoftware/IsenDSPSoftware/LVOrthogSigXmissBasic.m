function LVOrthogSigXmissBasic(N,F,A,B)
% LVOrthogSigXmissBasic(20,7,2,5)
% N is carrier cycle length in samples, F is carrier frequency, A and B
% are real numbers
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
n = 0:1:N-1;
C1 = cos(2*pi*n*F/N);
C2 = sin(2*pi*n*F/N);
S = A*C1 - B*C2;
Sr1 = (2/N)*sum(S.*C1)
Sr2 = -(2/N)*sum(S.*C2)



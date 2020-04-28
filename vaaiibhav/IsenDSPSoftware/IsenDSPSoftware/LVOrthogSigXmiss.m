function [Sr1,Sr2] = LVOrthogSigXmiss(N,F,A,B)
% N is carrier cycle length in samples, F is carrier frequency, A and B
% are equal length row vectors of real numbers
% [Sr1,Sr2] = LVOrthogSigXmiss(20,7,[2,-1,3,0,7],[5,2,-6,3,1])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(length(A)==length(B))
error('A and B must be the same length'); end
n = 0:1:N-1;
C1 = cos(2*pi*n*F/N);
C2 = sin(2*pi*n*F/N);
C1Mat = (C1')*(ones(1,length(A)));
C2Mat = (C2')*(ones(1,length(B)));
AMat = ones(N,1)*A; BMat = ones(N,1)*B;     
S1 = C1Mat.*(AMat); S1 = S1(:);
S2 = C2Mat.*(BMat); S2 = S2(:);
S = S1 - S2;
% must break S into one cycle frames
SigMat = reshape(S,N,fix(length(S)/N));
Sr1 = (2/N)*sum(SigMat.*C1Mat);
Sr2 = -(2/N)*sum(SigMat.*C2Mat);

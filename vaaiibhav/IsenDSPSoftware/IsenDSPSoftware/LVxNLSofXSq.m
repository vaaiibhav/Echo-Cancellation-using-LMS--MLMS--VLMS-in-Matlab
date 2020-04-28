function [yC,nC] = LVxNLSofXSq(NLCoeff,x)
% Delays, weights, and sums the square of the
% input sequence x ( = x[n] with n = 0:1:length(x)-1)
% according to yC = c(0)*x(n).^2 + c(1)*x(n-1).^2 + ... 
% with NLCoeff (= [c[0],c[1],c[2],...]) and
% nC are the sample indices of yC.
% Test call:
% [yC,nC] =  LVxNLSofXSq([1,-2,1],cos(2*pi*12*[0:1:63]/64) )
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
x1 = NLCoeff(1)*(x.^2); 
nC = [0:1:length(x1)-1];
yC = x1;
if length(NLCoeff)<2
    return
end
for NLCoeffCtr = 2:1:length(NLCoeff)
xC = [NLCoeff(NLCoeffCtr)*(x.^2)];
newnC = [0:1:length(xC)-1] + (NLCoeffCtr-1);
[yC, nC] = LVAddSeqs(yC,nC,xC,newnC);
end
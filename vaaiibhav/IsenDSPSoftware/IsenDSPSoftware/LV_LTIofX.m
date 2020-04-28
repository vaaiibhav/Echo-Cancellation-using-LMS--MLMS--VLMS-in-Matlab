function [yC,nC] = LV_LTIofX(LTICoeff,x)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
% [yC,nC] = LV_LTIofX([1,-2,1],cos(2*pi*12*[0:1:63]/64) )
x1 = LTICoeff(1)*x; 
nC = [0:1:length(x1)-1];
yC = x1;
if length(LTICoeff)<2
    return
end
for LTICoeffCtr = 2:1:length(LTICoeff)
xC = [LTICoeff(LTICoeffCtr)*x];
newnC = [0:1:length(xC)-1] + (LTICoeffCtr-1);
[yC, nC] = LVAddSeqs(yC,nC,xC,newnC);
end
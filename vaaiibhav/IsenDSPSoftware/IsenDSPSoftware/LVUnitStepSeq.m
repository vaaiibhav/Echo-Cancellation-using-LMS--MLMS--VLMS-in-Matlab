function [yVals,xIndices] = LVUnitStepSeq(n,Nlow,Nhigh)
% [yVals,xIndices] = LVUnitStepSeq(-2,-10,10)
% n is location of the beginning of the unit step in the interval of samples indices running
% from Nlow to Nhigh, i.e., Nlow:1:Nhigh
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
xIndices = [Nlow:1:Nhigh];
yVals(1:1:length(xIndices)) = 0;
posZInd = find((xIndices-n)==0);
yVals(posZInd:1:length(yVals)) = 1; 
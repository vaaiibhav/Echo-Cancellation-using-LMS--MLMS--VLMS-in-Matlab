function [xVals,xIndices] = LVUnitImpSeq(n,Nlow,Nhigh)
% [xVals,xIndices] = LVUnitImpSeq(-2,-10,10)
% n is location of unit impulse in the interval of samples indices running
% from Nlow to Nhigh, i.e., Nlow:1:Nhigh
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

xIndices = [Nlow:1:Nhigh];
xVals = [(xIndices-n) == 0];
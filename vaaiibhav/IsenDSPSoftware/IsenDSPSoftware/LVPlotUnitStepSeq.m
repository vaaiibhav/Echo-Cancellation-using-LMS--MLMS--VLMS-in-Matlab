function LVPlotUnitStepSeq(n,Nlow,Nhigh)
% LVPlotUnitStepSeq(-2,-20,100)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    xIndices = [Nlow:1:Nhigh];
    yVals(1:1:length(xIndices)) = 0;
    posZInd = find((xIndices-n)==0)
    yVals(posZInd:1:length(xIndices)) = 1; 
    stem(xIndices,yVals)
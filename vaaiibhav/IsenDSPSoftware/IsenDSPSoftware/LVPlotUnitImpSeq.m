function LVPlotUnitImpSeq(n,Nlow,Nhigh)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    xIndices = [Nlow:1:Nhigh]; 
    xVals = zeros(1,length(xIndices));
    xVals(find(xIndices-n==0))=1;
    stem(xIndices,xVals)
    

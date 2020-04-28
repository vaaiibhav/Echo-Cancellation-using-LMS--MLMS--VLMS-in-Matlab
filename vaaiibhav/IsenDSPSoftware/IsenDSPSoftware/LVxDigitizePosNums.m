function [OutputMat] = LVxDigitizePosNums(xNumsToDig,NumBits)
% function [OutputMat] = LVxDigitizePosNums(xNumsToDig,NumBits)
% Receives a vector of positive decimal numbers to convert to 
% binary notation using NumBits number of bits, and produces as
% Output a matrix each row of which is a binary representation of the
% corresponding input decimal number.
% Test call:
% LVxDigitizePosNums([0,1,34,23,2,17,254,255,127],8)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
xWtVec = 2.^(NumBits-1:-1:0);
xhowManySamps = length(xNumsToDig);
OutputMat = zeros(xhowManySamps,NumBits);
for BitCtr = 1:1:NumBits 
    OutputMat(1:xhowManySamps,BitCtr) = 1; 
    TestMat = OutputMat*xWtVec' - xNumsToDig';   
    [i,j] = find(TestMat>0);
    OutputMat(i,BitCtr) = 0;
 end
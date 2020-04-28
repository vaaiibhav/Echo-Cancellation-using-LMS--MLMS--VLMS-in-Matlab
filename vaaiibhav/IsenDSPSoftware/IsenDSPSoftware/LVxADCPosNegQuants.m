function  BinaryOutput = LVxADCPosNegQuants(NumsToDig)
%function  BinaryOutput = LVxADCPosNegQuants(NumsToDig)
%Returns a matrix of binary numbers in sign plus
%magnitude format constituting the binary conversion of the decimal numbers
%input as the vector NumsToDig
% Test call:
% [BinOut] = LVxADCPosNegQuants([17 7 -100  77  -63  98  21  -35])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
NumsToDig = real(NumsToDig);
LargestMagnitude = max(abs(NumsToDig));      
n = 1;
while ((2^n)-1)<LargestMagnitude
   n = n+1;
end
NumberOfBits = n + 1; 
[nRows,nCols] = size(NumsToDig);
if nRows==1 & nCols>1
   howManySamps = nCols; 
elseif nRows > 1 & nCols==1
   NumsToDig = NumsToDig';
   howManySamps = nRows;
elseif nRows==1 & nCols==1
   howManySamps=1;
elseif nRows>1 & nCols>1
	Comment = 'Improper format; ending procedure'
   return
end
NegIndices = find(NumsToDig<0);
NumsToDig = abs(NumsToDig);

WtVec(1,1:NumberOfBits-1) = 2.^(0:1:NumberOfBits-2); % vectorized 

OutputLessSignBit = DigitizePosNums(NumsToDig,NumberOfBits-1);
BinaryOutput(1:howManySamps,1) = zeros(howManySamps,1);
BinaryOutput(1:howManySamps,2:NumberOfBits) = OutputLessSignBit;
BinaryOutput(NegIndices,1) = ones(length(NegIndices),1);

function [OutputMat] = DigitizePosNums(xNumsToDig,NumBits) 
xWtVec = 2.^(NumBits-1:-1:0);
xhowManySamps = length(xNumsToDig);
OutputMat = zeros(xhowManySamps,NumBits);
for BitCtr = 1:1:NumBits 
    OutputMat(1:xhowManySamps,BitCtr) = 1; 
    TestMat = OutputMat*xWtVec' - xNumsToDig';   
    [i,j] = find(TestMat>0);
    OutputMat(i,BitCtr) = 0;
 end

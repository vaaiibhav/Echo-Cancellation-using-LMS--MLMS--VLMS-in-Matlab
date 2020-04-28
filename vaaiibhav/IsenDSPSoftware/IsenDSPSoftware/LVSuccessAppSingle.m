function [BinOut, Err] = LVSuccessAppSingle(DecNum,Bits) 
% [BinMat, Err] = LVSuccessAppSingle(5.75,4)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if length(DecNum)>1; error('DecNum must be scalar')
end
BinOut = zeros(1,Bits); BinWord = zeros(1,Bits);
for ctr = 1:1:Bits
Status = (['Setting Bit ',num2str(ctr),' to 1'])
BinWord(1,ctr) = 1
DecEquivCurrBinWord = LVBinary2DecimalVec(BinWord)
Status = (['Subtracting decimal equiv of BinWord from sample being quantized'])
diff = DecNum - DecEquivCurrBinWord
if diff < 0;
Status = (['Resetting Bit ',num2str(ctr),' to 0'])
BinWord(1,ctr) = 0; 
end
end
Status = (['Final Binary Word:'])
BinOut(1,:) = BinWord;
Err = DecNum - LVBinary2DecimalVec(BinWord);
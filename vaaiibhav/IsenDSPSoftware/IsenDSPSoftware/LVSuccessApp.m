function [BinMat,Err] = LVSuccessApp(DecNum,MaxBits,LSBBias) 
% [BinMat,Err] = LVSuccessApp(7.5*[sin(0.5*pi*[0:1/18:1])]+7.5,4,1);
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    BinMat = zeros(length(DecNum),MaxBits);
    DecOutMat = zeros(1,length(DecNum));
    if LSBBias==0; xDecNum = DecNum;
    else; xDecNum = DecNum + 0.5; end
    for DecNumCtr = 1:1:length(DecNum)
    BinTstWord = zeros(1,MaxBits);
     for ctr = 1:1:MaxBits
     BinTstWord(1,ctr) = 1;
     DecEquivCurrBinTstWord = LVBinary2DecimalVec(BinTstWord);
     diff = xDecNum(DecNumCtr) - DecEquivCurrBinTstWord;
     if diff < 0; BinTstWord(1,ctr) = 0;
     elseif diff==0; break; else; end
     end 
     BinMat(DecNumCtr,:) = BinTstWord;
     DecOutMat(1,DecNumCtr) = LVBinary2DecimalVec(BinTstWord);
     Err = DecNum-DecOutMat;
    end
    figure(78); clf; subplot(211); ldn = length(DecNum);
    xvec = 0:1:ldn-1; hold on; plot(xvec,DecNum,'b:'); 
    plot(xvec,DecNum,'bo'); stairs(xvec,DecOutMat); 
    xlabel('(a) Sample'); ylabel('Amplitude')
    subplot(212); stairs(xvec,Err); 
    xlabel('(b) Sample'); ylabel('Error')
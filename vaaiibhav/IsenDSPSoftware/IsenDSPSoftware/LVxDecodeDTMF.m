function [DialedDigits] = LVxDecodeDTMF(InputSignal)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

% Divide signal into windows of length 0.25 second, separated by 0.125
% second

NumDigits = round(length(InputSignal)/3000);

%RowMat = [697, 240, 21;  770, 374, 36;  852, 338, 36; 941,170,20];
%ColMat = [1209, 258, 39 ; 1336, 461, 77 ; 1477, 222, 41  ; 1633, 191, 39 ];

RowMat = [697, 200, 17;  770, 200, 19;  852, 200, 21; 941,200,24];
ColMat = [1209, 200, 30 ; 1336, 200, 33 ; 1477, 200, 37 ; 1633, 200, 41 ];

for ctr = 1:1:NumDigits
    
WinStart = 3000*(ctr-1)+1;
WinEnd = 3000*ctr-1000;

TestSig = InputSignal(WinStart:WinEnd);

for rowctr = 1:1:4
    testFreq = RowMat(rowctr,1);
    WinLen = RowMat(rowctr,2);
    Bin = RowMat(rowctr,3);
    netTestSig = InputSignal(WinStart+1+500:WinStart+WinLen+500);
    [DFTGBin,MagDFTGBin] = LVxDFTViaGoertzelBin(Bin,netTestSig);
    RowMag(rowctr) = (10/WinLen)*MagDFTGBin;
end

[theLocs, thePeaks] = LVfindPeaks(RowMag,1,1);
theRow = theLocs(1);

for colctr = 1:1:4
    testFreq = ColMat(colctr,1);
    WinLen = ColMat(colctr,2);
    Bin = ColMat(colctr,3);
    netTestSig = InputSignal(WinStart+1:WinStart+WinLen);
    [DFTGBin,MagDFTGBin] = LVxDFTViaGoertzelBin(Bin,netTestSig); 
    ColMag(colctr) = MagDFTGBin;
end
[theLocs, thePeaks] = LVfindPeaks(ColMag,1,1);
theCol = theLocs(1);

DigitIndices(ctr,1) = theRow;
DigitIndices(ctr,2) = theCol;
end

% ans = DigitIndices
OutputDigits = [];
a = size(DigitIndices);
for ctr = 1:1:a(1)
if DigitIndices(ctr,1)==1
    if DigitIndices(ctr,2)==1
        OutputDigits = [OutputDigits,' ',num2str(1)];
    elseif DigitIndices(ctr,2)==2
        OutputDigits = [OutputDigits,' ',num2str(2)];
    elseif DigitIndices(ctr,2)==3
        OutputDigits = [OutputDigits,' ',num2str(3)];
    elseif DigitIndices(ctr,2)==4
        OutputDigits = [OutputDigits,' ','A'];
    end
elseif DigitIndices(ctr,1)==2
      if DigitIndices(ctr,2)==1
        OutputDigits = [OutputDigits,' ',num2str(4)];
    elseif DigitIndices(ctr,2)==2
        OutputDigits = [OutputDigits,' ',num2str(5)];
    elseif DigitIndices(ctr,2)==3
        OutputDigits = [OutputDigits,' ',num2str(6)];
    elseif DigitIndices(ctr,2)==4
        OutputDigits = [OutputDigits,' ','B'];
    end
elseif DigitIndices(ctr,1)==3
          if DigitIndices(ctr,2)==1
        OutputDigits = [OutputDigits,' ',num2str(7)];
    elseif DigitIndices(ctr,2)==2
        OutputDigits = [OutputDigits,' ',num2str(8)];
    elseif DigitIndices(ctr,2)==3
        OutputDigits = [OutputDigits,' ',num2str(9)];
    elseif DigitIndices(ctr,2)==4
        OutputDigits = [OutputDigits,' ','C'];
    end
elseif DigitIndices(ctr,1)==4
          if DigitIndices(ctr,2)==1
        OutputDigits = [OutputDigits,' ','*'];
    elseif DigitIndices(ctr,2)==2
        OutputDigits = [OutputDigits,' ',num2str(0)];
    elseif DigitIndices(ctr,2)==3
        OutputDigits = [OutputDigits,' ','#'];
    elseif DigitIndices(ctr,2)==4
        OutputDigits = [OutputDigits,' ','D'];
    end
end
end
DialedDigits = OutputDigits;
% OW = EncodeDTMF([7 0 3 ],0.05)
% [DialedDigits] = DecodeDTMF(OW)
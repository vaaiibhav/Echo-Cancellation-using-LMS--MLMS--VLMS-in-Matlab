function OutMat = LVxVector2FramesInMatrix(Sig,SzWin,SampsOvrLap)
% function OutMat = LVxVector2FramesInMatrix(Sig,SzWin,SampsOvrLap)
% Divides an input signal vector Sig into frames of length SzWin, with 
% an amount of overlap in samples equal to SampsOvrLap. The output matrix
% Each column of OutMat is one frame of the input Sig.
% Test calls:
% OutMat = LVxVector2FramesInMatrix([0:1:33],8,4)
% OutMat = LVxVector2FramesInMatrix([0:1:33],8,0)
% OutMat = LVxVector2FramesInMatrix([0:1:33],8,1)
% OutMat = LVxVector2FramesInMatrix([0:1:33],8,2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
lenFile = length(Sig)
if SampsOvrLap==0
nMinus1 = ceil(lenFile/SzWin)-1;
NoFrames = nMinus1 + 1;
totSampsNeeded = SzWin*NoFrames;
else
nMinus1 = ceil((lenFile-SzWin)/SampsOvrLap);
NoFrames = nMinus1 + 1;
totSampsNeeded = SzWin + NoFrames*SampsOvrLap;
end

Sig = Sig(:)';
Sig = [Sig,zeros(1,totSampsNeeded-lenFile)];
Sig = Sig';

if SampsOvrLap==0
for n = 1:1:NoFrames
OutMat(1:SzWin,n) = Sig( 1+(n-1)*SzWin : SzWin +(n-1)*SzWin );
end    
    
else
    
for n = 1:1:NoFrames
OutMat(1:SzWin,n) = Sig( 1+(n-1)*SampsOvrLap:SzWin+(n-1)*SampsOvrLap );
end

end

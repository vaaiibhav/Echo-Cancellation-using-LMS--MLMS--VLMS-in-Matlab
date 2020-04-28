function tdSig = LVxTDMat2SigVec(TDMat,SzWin,SampsOverLap)
% Generates a signal vector tdSig from a matrix of time domain frames of a 
% signal, each frame of length SzWin, with SampsOverLap samples of overlap. 
% The code determines which dimension of TDMat matches SzWin and assumes 
% that the other dimension is the frame index. If neither dimension of 
% TDMat matches SzWin, an error is thrown. 
%
% Test call pair for OverLap=0
% OutMat = LVxVector2FramesInMatrix([0:1:33],8,0)
% tdSig = LVxTDMat2SigVec(OutMat,8,0)
%
% The following four lines of code, when run in sequence, return the input vector [0:1:33]:
% OutMat = LVxVector2FramesInMatrix([0:1:33],8,4); 
% szO = size(OutMat);
% OutMat = OutMat.*(triang(8)*ones(1,szO(2)))
% tdSig = LVxTDMat2SigVec(OutMat,8,4)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
szTDMat = size(TDMat);
if szTDMat(1)== SzWin
   IntNoWin = szTDMat(2);
   TDMat = TDMat';  % orient so each row is a frame for code below
elseif szTDMat(2)==SzWin
    IntNoWin = szTDMat(1);
else
    error('One dimension of TDMat must match SzWin')  
end

tdSig = [];
if SampsOverLap==0
lenFile = IntNoWin*SzWin;    
else
lenFile = 1 + SampsOverLap*(IntNoWin) + SzWin;
end
tdSig(1,1:lenFile + 2*SzWin) = 0;
sztdSig = size(tdSig);

tdSig(1,1:SzWin)= tdSig(1,1:SzWin)+ TDMat(1,1:SzWin);

for n = 1:1:IntNoWin-1;
   tdSig(1,(n*SzWin + 1 -(n)*SampsOverLap):((n+1)*SzWin-(n)*SampsOverLap))...
  = tdSig(1,(n*SzWin + 1 -(n)*SampsOverLap):((n+1)*SzWin-(n)*SampsOverLap)) ...
      + TDMat(n+1,1:SzWin);
end
tdSig = tdSig';
tdSig = tdSig(:)';
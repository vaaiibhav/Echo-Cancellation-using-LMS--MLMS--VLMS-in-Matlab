function LVManualLPFRemezExch(L,LenGrid,wp,ws,curXFrqs)
% LVManualLPFRemezExch(9,145,0.45,0.55,[0,0.225,0.45,0.55,0.775,1])
% LVManualLPFRemezExch(19,304,0.45,0.55,[0,0.142,0.27,0.385,0.45,0.55,0.59,0.68,0.78,0.89,1])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
     NormFrGrid = [0:1:LenGrid-1]/(LenGrid-1);
     FrGrid = pi*NormFrGrid; WtVec = ones(1,LenGrid);
     wc = (wp+ws)/2; XFrindOnFG = round(curXFrqs*(LenGrid-1)+1);
     Hdr = ones(1,round(wc*(LenGrid-1)));
     Hdr = [Hdr, zeros(1,LenGrid-length(Hdr))];
     Q = 1; kLim = (L-1)/2; WMat(1:kLim+2,1) = 1; Hdr = Hdr./Q;
     WtVec = WtVec.*Q; k = 0:1:kLim+1;
     Num = ((-1).^k); Denom = WtVec(1,XFrindOnFG);
     WMat(k+1,2:kLim+1) = cos(pi*curXFrqs(1,k+1)'*([1:kLim]));
     WMat(k+1,kLim+2) = (Num./Denom)';
     HdrVec = Hdr(round((LenGrid-1)*[curXFrqs]+1));
     AlDelVec = pinv(WMat)*(HdrVec)';
     delta = AlDelVec(length(AlDelVec)),
     P = 0; for pCtr = 1:length(AlDelVec)-1;
     P = P + AlDelVec(pCtr)*(cos(FrGrid*(pCtr-1)) ); end;
     E = WtVec.*([Hdr - P]); figure(8); clf;
     subplot(211); hold on; ad = abs(delta);
     loBndX = [0:1:round((LenGrid-1)*wp)]/(LenGrid-1);
     hiBndX = [round((LenGrid-1)*ws):1:LenGrid-1]/(LenGrid-1);
     plot(loBndX,E(1,1:round(wp*(LenGrid-1))+1));
     line([0,max(loBndX)],[ad, ad]); line([0,max(loBndX)],[-ad, -ad]);
     xlabel('Norm Freq'); ylabel('Err (Passband)');
     subplot(212); hold on;
     plot(hiBndX,E(1,round(ws*(LenGrid-1))+1:LenGrid));
     line([ws,1],[ad, ad]); line([ws,1],[-ad, -ad]);
     xlabel('Norm Freq'); ylabel('Err (Stopband)');
     Hr = Q.*P; figure(9); LenHr = length(Hr);
     plot([0:1:LenHr-1]/LenHr,Hr);
     xlabel('Norm Freq'); ylabel('Amp')
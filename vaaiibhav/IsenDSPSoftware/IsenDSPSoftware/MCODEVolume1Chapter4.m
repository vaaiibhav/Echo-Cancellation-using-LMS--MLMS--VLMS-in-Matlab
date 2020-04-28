% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% this script is not meant to be called; it contains clean copies of m-code
% from the text. To run any of the example m-code segments below, copy the
% code, paste it into the Command line, and press Return.
return
% Page 126-----------------------------------------------------------------
[CorC] = LVCorrCosinesZerothLag(1,17,16)

[CorC] = LVCorrCosinesZerothLag(0,16,16)

[CorC] = LVCorrCosinesZerothLag(8,40,16)

[CorC] = LVCorrCosinesZerothLag(6,55,16)
% Page 128----------------------------------------------------------------
s = LVSumCE(2,32)
s = LVSumCE(32,32)
% Page 130-----------------------------------------------------------------
LVxTestReconSineVariablePhase(1,32,-145)
% Page 131-----------------------------------------------------------------
k = 0;C=sum(cos(2*pi*k*(0:1:7)/8).^2)
% Page 132-----------------------------------------------------------------
LVxFreqTest(5,32,[],1)
% Page 136-----------------------------------------------------------------
LVOrthogSigXmissBasic(20,7,2,5)
% -------------------------------------------------------------------------
[Sr1,Sr2] = LVOrthogSigXmiss(20,7,[2,-1,3,0,7],[5,2,-6,3,1])
%Page 138------------------------------------------------------------------
ML_Correlation
%Page 139------------------------------------------------------------------
y = xcorr([sin(2*pi*(0:1:7)/8)],[cos(2*pi*(0:1:7)/8)])
figure; stem(y)
% -------------------------------------------------------------------------
y = xcorr([sin(2*pi*2*(0:1:7)/8)]); figure; stem([-7:1:7],y)
%Page ---------------------------------------------------------------------
LVCorrSeqSinOrthog(32,128,1,8,0)
%Page 140------------------------------------------------------------------
LVCorrSeqSinOrthog(32,128,1,12,0)
%Page 143------------------------------------------------------------------
ycorr = xcorr([5,4,3,2,1],[1,2,3,4,5])
%Page 143------------------------------------------------------------------
yconv = conv([5,4,3,2,1], fliplr([1,2,3,4,5]))
%Page 145------------------------------------------------------------------
LVxMatchedFilter(0.5,128,1)
%Page 146------------------------------------------------------------------
LVFreqResp([1.9,0,-0.9,0,1.9,0,-3.2],500)
%Page 147------------------------------------------------------------------
LVFreqResp([1.9, 0, -0.9, 0, 1.9, 0, -3.2], 500)
%Page 150------------------------------------------------------------------
N=54; k = 9; x = cos( 2*pi*k*( 0:1:N-1 )/N);
LVFreqResp(x, 500)
%Page 150------------------------------------------------------------------
Imp = LVBasicFiltMultCorr(31,0,7)
Imp = LVBasicFiltMultCorr(30,0,4);
Imp = LVBasicFiltMultCorr(30,5,9);
Imp = LVBasicFiltMultCorr(30,10,15);
%Page 155------------------------------------------------------------------
LVxCorrDelayMeasure(0.1)
%Page 159------------------------------------------------------------------
SR = 24; b = 1; p = 0.8; y = zeros(1,SR);
x = [1,zeros(1,SR)]; y(1) = b*x(1);
for n = 2:1:SR
y(1,n) = b*x(1,n) + p*y(1,n - 1);
end; figure;
stem(y)
%Page 161------------------------------------------------------------------
p = 0.9*j; y = filter([1], [1,-p], [ones(1,150)])
% -------------------------------------------------------------------------
p = 0.9*j; xp = 0:1:99; x = p.^xp; y = conv(x, ones(1,200))
%Page 164------------------------------------------------------------------
x = ones(1,100); y = filter(0.1,[1,-0.9],x); figure;stem(y)
% ------------------------------------------------------------------------
x = 0.1*ones(1,100); y = filter(1,[1,-0.9],x); figure;stem(y)
%Page 166------------------------------------------------------------------
ML_SinglePole
%Page 169------------------------------------------------------------------
LVImpCmpxConjPoles(0.9,24)
%Page 170------------------------------------------------------------------
x = [1,zeros(1,50)]; y = filter(1,[1,-1.8,0.81],x);
figure; stem(y)
%Page 171------------------------------------------------------------------
y = conv([1, -(0.672 +j*0.672)],[1, -(0.672 - j*0.672)])
% -------------------------------------------------------------------------
LVRealFiltfromCCPoles(0.95, pi/4)
% -------------------------------------------------------------------------
LVRealFiltfromCCPoles(0.95,pi/4)
% -------------------------------------------------------------------------
ML_DragPoleZero
% -------------------------------------------------------------------------






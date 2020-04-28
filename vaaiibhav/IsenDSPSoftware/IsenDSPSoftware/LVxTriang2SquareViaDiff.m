function LVxTriang2SquareViaDiff(SR,HiHarmNo,TrFundFreq,Ldiff)
% LVxTriang2SquareViaDiff(2000,81,20,94)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
n = 0:1:SR; kVec = 1:2:HiHarmNo;
AmpVec = (1./(kVec.^2))';
Cmat = ((2*pi*n*TrFundFreq/SR)')*kVec;
y = cos(Cmat)*AmpVec;
Imp = LVxDifferentiatorTypeIV(Ldiff);
derivTriang =conv(Imp,y); 
del = fix( length(Imp)/2);
figure(8); subplot(211); plot(y(20:400)); 
xlabel('Triangle Waveform'); subplot(212); 
plot(derivTriang(20+del:400+del))
xlabel('Derivative of Triangle Waveform')
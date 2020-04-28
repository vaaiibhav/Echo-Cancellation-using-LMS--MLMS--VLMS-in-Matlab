function LVCombFilter(Tau)
% LVCombFilter(5)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
ImpAdd = [1 zeros(1,Tau-1) 1];
ImpSub = [1 zeros(1,Tau-1) -1];
DTFTLen = 1024; xplot = [1:1:DTFTLen/2+1];
xvec= (xplot-1)/(DTFTLen/2);
subplot(2,1,1); yAdd = abs(fft(ImpAdd,DTFTLen));
plot(xvec,yAdd(xplot),'b'); ylabel(['Magnitude'])
xlabel(['(a)  Normalized Frequency'])
axis([0 1 0 1.2*max(yAdd)])
subplot(2,1,2); ySub = abs(fft(ImpSub,DTFTLen));
plot(xvec,ySub(xplot),'b'); ylabel(['Magnitude'])
xlabel(['(b)  Normalized Frequency'])
axis([0  1 0 1.2*max(ySub)])



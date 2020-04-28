function LVxFreqSampFilter(Imp,CFsB,CFsA,BFs,AFs,x)
% LVxFreqSampFilter(Imp,CFsB,CFsA,BFs,AFs,x)
% Receives an impulse response Imp and a signal x, filters x and displays
% the result two different ways, first, using the impulse response itself
% as a Direct Form FIR, and second, using Frequency Sampling implementation
% coefficients CFsB,CFsA,BFs,AFs. 
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
N = length(Imp);
yFs = filter(CFsB,CFsA,x);

szB = size(BFs);
NoFiltBsec = szB(1);
szA = size(AFs);
NoFiltAsec = szA(1);

if ~(szA==szB)
    error('Number of filter sections defined by first dimension of Bfs must be same as defiend by first dimension of Afs')
end
y = 0;
for ctr = 1:1:NoFiltBsec   
    y =  y + filter(BFs(ctr,:),AFs(ctr,:),yFs);
end

FsFilt = y;

ImpFilt = filter(Imp,1,x);

figure(88)

subplot(211)

plot(ImpFilt)
xlabel('(a) Sample')
ylabel('Amplitude')


subplot(212)

plot(FsFilt)

xlabel('(b) Sample')
ylabel('Amplitude')



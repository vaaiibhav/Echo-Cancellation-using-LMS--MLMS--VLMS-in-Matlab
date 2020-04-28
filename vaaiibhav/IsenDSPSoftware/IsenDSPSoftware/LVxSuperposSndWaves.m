function LVxSuperposSndWaves(minFreq,maxFreq,FreqInc,d,A)
% LVxSuperposSndWaves(0,1000,10,0.001,0.9)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
fVec = (minFreq:FreqInc:maxFreq)';
t = 0.99:0.0005:1.01; 
f = zeros(length(fVec),2); 
f(:,1) = fVec;

for ctr = 1:1:length(fVec)
y = exp(j*2*pi*t*fVec(ctr)) + A*exp(j*2*pi*(t-d)*fVec(ctr));
f(ctr,2) = max(abs(y));
end

ffs = f(:,2);
plotlim = 1.2*max(abs(ffs));

figure(4)
plot(f(:,1),f(:,2))

axis([minFreq,maxFreq,0,plotlim])

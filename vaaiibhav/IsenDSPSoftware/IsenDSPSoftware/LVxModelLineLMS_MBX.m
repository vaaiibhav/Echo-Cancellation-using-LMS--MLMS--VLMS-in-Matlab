function LVxModelLineLMS_MBX(M, B, xVec, Mu, bMu2mMuRat,NoIts)
% function LVxModelLineLMS_MBX(M, B, xVec, Mu, bMu2mMuRat, NoIts)
%
% xVec and yVec specify a set of points in a plane having the corresponding 
% x and y coordinates to be modeled by the line y = mx + b where m is the
% slope and b is the y-intercept. The parameters m and b are determined
% using an adaptive process having an update weighting value of Mu which is
% conducted for NoIts iterations.
%
% Typical call:
%
% LVxModelLineLMS_MBX(3,-1,[0 1 2 3 ],1.2, 3,30)
% LVxModelLineLMS_MBX(3,-1,[-10:1:10],0.5, 30,5)
% LVxModelLineLMS_MBX(3,-1,[0:0.1:2],1.15,3,45)
% LVxModelLineLMS_MBX(3,-1,[0:1:2],1.1,2,30)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

N = length(xVec);
m = zeros(1,NoIts);
b = zeros(1,NoIts);
MSE = zeros(1,NoIts);
plotxVec = min(xVec):0.1:max(xVec);

yVec = M*xVec + B;

Mu = Mu/max(abs(yVec));

figure(888)
clf

for Ctr = 1:1:NoIts
MSE(Ctr) = (1/N)*sum( (yVec -(m(Ctr)*xVec + b(Ctr))).^2 );
PartDerivM = (2/N)*sum( (yVec -(m(Ctr)*xVec + b(Ctr))).*(-xVec));
PartDerivB = (2/N)*sum( (yVec -(m(Ctr)*xVec + b(Ctr)))*(-1));
m(Ctr+1) = m(Ctr) - Mu*PartDerivM;
b(Ctr+1) = b(Ctr) - bMu2mMuRat*Mu*PartDerivB;
 
if Ctr>1
    if MSE(Ctr)>1.02*MSE(Ctr-1)
        m(Ctr+1) = m(Ctr);
        b(Ctr+1) = b(Ctr);
        Mu = Mu/2
    end
end
end

finalB = b(NoIts)
finalM = m(NoIts)

plotyVec = m(Ctr+1)*plotxVec + b(Ctr+1);

xRange = abs(max(xVec) - min(xVec));
maxx = max(xVec) + 0.05*xRange;
minx = min(xVec) - 0.05*xRange;

miny = min( [min(yVec)  min(plotyVec)] );
maxy = max( [max(yVec)  max(plotyVec)] );

yRange = abs(maxy  -  miny);
Netmaxy = maxy + 0.1*yRange;
Netminy = miny - 0.1*yRange;

subplot(211)

plot(xVec,yVec,'b.')
hold on
plot(plotxVec,plotyVec)
grid on
xlabel(['(a) x'])
ylabel('y')
axis([  minx  maxx  Netminy  Netmaxy ])
title(['Final coefficient estimates: m = ',num2str(finalM),'; b = ',num2str(finalB)])

subplot(212)

plot(20*log10(MSE))
grid on
xlabel('(b) Iteration')
ylabel('MSE, dB')
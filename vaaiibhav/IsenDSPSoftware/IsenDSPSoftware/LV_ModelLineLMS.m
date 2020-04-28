function LV_ModelLineLMS(xVec,yVec,Mu, NoIts)
% function LV_ModelLineLMS(xVec,yVec,Mu, NoIts)
% Receives a set of test points in a plane as 
% xVec and yVec, a value of Mu suitable for use 
% in a gradient-based coefficient update algorithm,
% and a number of iterations NoIts to perform
% in an attempt to produce values of slope and 
% y-intercept that best model the data as a straight
% line.
% Test call:
% LV_ModelLineLMS([0,1,2,3],[-1,1,3,5],0.05,50)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(length(xVec) == length(yVec))
    Comment = 'xVec and yVec must have same length!'
    return
end

N = length(xVec);
m = zeros(1,NoIts);
b = zeros(1,NoIts);
MSE = zeros(1,NoIts);
plotxVec = min(xVec):0.1:max(xVec);

scnsize = get(0,'ScreenSize');
pos1 = [0.005*scnsize(3), 0.05*scnsize(4), 0.99*scnsize(3), 0.84*scnsize(4)];

figure(888)
clf
set(888,'color',[1 1 1],'numbertitle','off','position',pos1,'name','Modeling data with a line')

for Ctr = 1:1:NoIts
MSE(Ctr) = (1/N)*sum( (yVec -(m(Ctr)*xVec + b(Ctr))).^2 );
PartDerivM = (2/N)*sum( (yVec -(m(Ctr)*xVec + b(Ctr))).*(-xVec));
PartDerivB = (2/N)*sum( (yVec -(m(Ctr)*xVec + b(Ctr)))*(-1));
m(Ctr+1) = m(Ctr) - Mu*PartDerivM;
b(Ctr+1) = b(Ctr) - Mu*PartDerivB;
hold on
plotyVec = m(Ctr+1)*plotxVec + b(Ctr+1);
plot(xVec,yVec,'bo')
plot(xVec,yVec,'b')
if Ctr==1|Ctr==2|Ctr==NoIts
plot(plotxVec,plotyVec,':')
end

set(gca,'fontsize',14)
xlabel(['x'],'fontsize',16)
ylabel('y','fontsize',16)

xRange = abs(max(xVec) - min(xVec));
maxx = max(xVec) + 0.05*xRange;
minx = min(xVec); % - 0.05*xRange;

miny = min( [min(yVec)  min(plotyVec)] );
maxy = max( [max(yVec)  max(plotyVec)] );

yRange = abs(maxy  -  miny);
Netmaxy = maxy + 0.1*yRange;
Netminy = miny - 0.1*yRange;

title(['Coeff estimates: m = ',num2str(m(Ctr+1)),'; b = ',num2str(b(Ctr+1)),' at Iter. ',num2str(Ctr)],'fontsize',16)
axis([  minx  maxx  Netminy  Netmaxy ])
pause(0.05)
end

finalB = b(NoIts)
finalM = m(NoIts)

plotyVec = m(Ctr+1)*plotxVec + b(Ctr+1);
plot(xVec,yVec,'b.')
hold on
plot(plotxVec,plotyVec,'b:')
set(gca,'fontsize',14)
xlabel(['x'],'fontsize',16)
ylabel('y','fontsize',16)
axis([  minx  maxx  Netminy  Netmaxy ])
title(['Final coefficient estimates: m = ',num2str(finalM),'; b = ',num2str(finalB)],'fontsize',16)






function LVxModelPtswLine(TestM,TestYInt,xTest,Mu0,yMu2MuRatio,MStart,yIntStart,NoIts)
%function LVxModelPtswLine(TestM,TestYInt,xTest,Mu0,yMu2MuRatio,MStart,yIntStart,NoIts)
%
% The line to be modeled is represented in the point-slope form, and thus the
% test line creation parameters are TestM (the slope), TestYInt (the test Y-intercept).
% and the vector xTest. Mu0 is the initial value used to weight the gradient estimate update term;
% yMu2MuRatio is a ratio by which to relatively weight the partial derivative of
% y-intercept relative to slope (see text); MStart is the initial estimate of slope to use;
% yIntStart is the initial estimate of y-intercept to use, and NoIts is the
% number of iterations to perform.  
% The script creates separate 3-D plots of the performance surface and the
% coefficient track, both plots having the same axis limits for all three
% axes, enabling comparison.
% 
% Sample calls:
% LVxModelPtswLine(2,0,[-10:1:10],0.014,37,-10,8,12) % fast convergence
% LVxModelPtswLine(2,0,[-10:1:10],0.005,1,-10,8,40) % sharp turn in gradient search
% LVxModelPtswLine(2,0,[-10:1:10],0.025,1,-10,8,60) % overshoot on slope
% LVxModelPtswLine(2,0,[-10:1:10],0.025,40,-10,8,30) % overshoot both slope & yInt
% LVxModelPtswLine(20,-90,[-10:1:10],0.005,10,-10,8,40)
% LVxModelPtswLine(2,0,[0:1:10],0.01,50,-10,8,80); 
% LVxModelPtswLine(1,-2,[0:1:10],0.01,50,-10,8,80);
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

figure(769)
 clf
set(769,'color',[1,1,1])
 
% subplot(121)
 
thisMSE = [];
psMSE = [];
bott = 10^(-8);
TestDelta = 0.005;
maxM = max([abs(MStart),abs(TestM)]);
maxyInt = max([abs(yIntStart),abs(TestYInt) ]);

Slope = -2*maxM:maxM/10:2*maxM; % used to create performance surface
yIntercept = -2*maxyInt:maxyInt/10:2*maxyInt;  

yTest = TestYInt + TestM*xTest;

x = xTest;
Lenx = length(x);
 psMSE = zeros(length(Slope),length(yIntercept));
yEst = zeros(1,length(x));

for testMCtr = 1:1:length(Slope)
for yIntCtr = 1:1:length(yIntercept) 
 yEst(1:Lenx) = yIntercept(yIntCtr) + Slope(testMCtr)*x;
 psMSE(testMCtr,yIntCtr) = GetMSE(yTest,yEst(1:Lenx));
end
end

surf(yIntercept,Slope,psMSE)
xlabel(['y-Intercept'])
ylabel(['Slope'])
zlabel(['MSE'])
axis( [ min(yIntercept),max(yIntercept),min(Slope),max(Slope),0,max(max(psMSE))])
%view([-95,30])

thisMSE = zeros(1,NoIts+1);
yIntNow = yIntStart;
SlopeNow = MStart;

yEst(1:Lenx) = yIntStart + MStart*x;
thisMSE(1) = GetMSE(yTest,yEst(1:Lenx)); 

Mu = Mu0;
DoubleTestDelta = 2*TestDelta;

for MCtr = 1:1:NoIts  % the algorithm
yIntMu = yMu2MuRatio*Mu;
% get partial deriv for yIntercept
yEstTrialPlus = (yIntNow(MCtr) + TestDelta) + SlopeNow(MCtr)*x;
TestMSEPlus = GetMSE(yTest,yEstTrialPlus);
yEstTrialMinus = (yIntNow(MCtr)-TestDelta) + SlopeNow(MCtr)*x;
TestMSEMinus = GetMSE(yTest,yEstTrialMinus);
PartialyIntercept = (TestMSEPlus - TestMSEMinus)/DoubleTestDelta;
% get partial deriv for Slope
yEstTrialPlus = yIntNow(MCtr) + (SlopeNow(MCtr) + TestDelta)*x;
TestMSEPlus = GetMSE(yTest,yEstTrialPlus);
yEstTrialMinus = yIntNow(MCtr) + (SlopeNow(MCtr) - TestDelta)*x;
TestMSEMinus = GetMSE(yTest,yEstTrialMinus);
PartialSlope = (TestMSEPlus - TestMSEMinus)/DoubleTestDelta;
%============================================================
yIntNow(MCtr+1) = yIntNow(MCtr) - PartialyIntercept*yIntMu;
SlopeNow(MCtr+1) = SlopeNow(MCtr) - PartialSlope*Mu;

yEst(1:Lenx) = yIntNow(MCtr+1) + SlopeNow(MCtr+1)*x;
thisMSE(MCtr+1) = GetMSE(yTest,yEst(1:Lenx));

if thisMSE(MCtr+1)>10^6
    Comment = ['MSE excessively high; ending procedure; MCtr = ',num2str(MCtr)]
   return
end

end

finalMSE = thisMSE(MCtr+1)
finEst_yInt = yIntNow(MCtr+1)
finEst_Slope = SlopeNow(MCtr+1)

figure(770)
% clf
 set(770,'color',[1,1,1])
%subplot(122)
%hold on
plot3(yIntNow,SlopeNow,thisMSE,'bo')
%line(yIntNow,SlopeNow,thisMSE,'color',[1,1,1])
grid on
xlabel(['y-Intercept'])
ylabel(['Slope'])
zlabel(['MSE'])
axis( [ min(yIntercept),max(yIntercept),min(Slope),max(Slope),0,max(max(thisMSE))]) 
 % ========================================================================
%view([-95,30])

function [OutputMSE] = GetMSE(MyyTest,yEst)
[OutputMSE] = (1/length(yEst))*sum((MyyTest - yEst).^2);
return


      

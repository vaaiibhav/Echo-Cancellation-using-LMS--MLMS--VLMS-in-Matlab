function LVxGradientViaCP(FirstValX,StepSize,NoiseAmp,UseWgtingFcn,deltaX,NoIts)
% function LVxGradientViaCP(FirstValX,StepSize,NoiseAmp,UseWeightingFcn,deltaX,NoIts)
% FirstValX is the initial value of x; StepSize determines how
% large the coefficient update term can be, i.e., how much x can change by
% per iteration; 
%
% NoiseAmp specifies the amplitude of any noise to be added (may be 0); 
%
% If UseWeightingFcn is 1, the update term is additionally
% weighted with the current value of the cost function (reduces effective
% stepsize near convergence); pass UseWeightingFcn as 0 if not wanted.
%
% deltaX is the small perturbation amount used to estimate
% the partial derivative;
%
% NoIts is the maximum number of iterations to perform.
%
% A typical call is
%
% LVxGradientViaCP(4,0.0294,0.01,1,0.01,20)
% LVxGradientViaCP(4,0.05,0.01,0,0.01,50)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool


figure(8767)
clf

x = zeros(1,NoIts);
x(1) = FirstValX;

if FirstValX>0
allX = -1.25*FirstValX-1:0.1:1.25*FirstValX+1;
else
allX = 1.25*FirstValX-1:0.1:-1.25*FirstValX+1;    
end

allY = allX.^2 + 1;
XdeltaX = 1;

for n = 1:1:NoIts
   x(n) = x(n) + NoiseAmp*(rand(1,1)-0.5);
   CostFcn(n) = (x(n)-deltaX/2)^2 + 1;
   
   if CostFcn(n)> abs(FirstValX*2500) %10^4
    ylabel(['Value of Cost Fcn'])
    xlabel([ 'At Iteration ',num2str(n),', X = ',num2str(x(n),5),', CostFcn = ',num2str(CostFcn(n)),': Ending!'])
    
 	if FirstValX>0
       aFirstValX = max([FirstValX  1]);
    	axis([-1.5*aFirstValX  1.5*aFirstValX  0 (1.5*(aFirstValX^2)+2)])
 	else
    	aFirstValX = min([FirstValX  -1]);
    	axis([1.5*aFirstValX  -1.5*aFirstValX  0 (1.5*(aFirstValX^2)+2)])
  	end

   Comment = 'MSE > 10,000!--Ending Procedure!'
   break
   end   

 hold on
 plot(x(n),CostFcn(n),'bo');
 plot(allX,allY,'b:');
 
 if FirstValX>0
   aFirstValX = max([FirstValX  1]);
 else
   aFirstValX = min([FirstValX  -1]);
end

if n==1
 text(-0.4*abs(allX(1)),1.4*(aFirstValX^2+1),['Press Any Key to Continue'])
end

 TestCostFcn = (((x(n) + deltaX/2))^2 + 1);
 NetCostFcn = (TestCostFcn + CostFcn(n))/2;
 PartialDerivCstFcn = (TestCostFcn -  CostFcn(n))/deltaX; 

if UseWgtingFcn==1
   x(n+1) = x(n) - PartialDerivCstFcn*StepSize*CostFcn(n);
else
   x(n+1) = x(n) - PartialDerivCstFcn*StepSize;
end
  
   ylabel(['Value of Cost Fcn'])
   xlabel([ 'Iter: ',num2str(n),', X = ',num2str(x(n),4),', CostFcn = ',num2str(NetCostFcn,4),'; StepSize = ',num2str(StepSize,4),'; Gradient = ',num2str(PartialDerivCstFcn)])

   if FirstValX>0
   aFirstValX = max([FirstValX  1]);
   axis([-1.5*aFirstValX  1.5*aFirstValX 0 (1.5*(aFirstValX^2)+2)])
 	else
   aFirstValX = min([FirstValX  -1]);
   axis([1.5*aFirstValX  -1.5*aFirstValX 0 (1.5*(aFirstValX^2)+2)])
   end

   pause
end


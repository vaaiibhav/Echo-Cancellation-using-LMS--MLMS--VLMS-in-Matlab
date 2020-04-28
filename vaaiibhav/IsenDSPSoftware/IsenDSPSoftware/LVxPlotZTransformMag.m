function LVxPlotZTransformMag(ZeroVec,PoleVec,rMin,rMax,NSamps)
% function LVxPlotZTransformMag(ZeroVec,PoleVec,rMin,rMax,NSamps)\n'])
%
% LVxPlotZTransformMag is intended for use with right-handed z-transforms,
% i.e., those whose ROC is the area outside of a circle having a radius just
% larger than the magnitude of the largest pole of the z-transform.
%
% ZeroVec is a row vector of zeros in the z-plane; if no zeros, pass as []
% PoleVec is a row vector of poles in the z-plane; if no poles, pass as []
% rMin and rMax are the radii of the smallest and largest contours to
% compute; rMin and rMax may be equal to compute one contour. Pass rMin = rMax = 1 
% for standard frequency response.
% NSamps is the number of frequencies at which to evaluate the z-transform along 
% each contour. rMin is set at 0.02 greater than
% the magnitude of the largest pole unless the pole magnitude is between
% 0.98 and 1.0 in which case rMin is set at 1.0.
%
% Sample calls: 
%
% LVxPlotZTransformMag([-1 1],[0.95*exp(j*2*pi/5) 0.95*exp(-j*2*pi/5)],1,1,256)
% LVxPlotZTransformMag([exp(j*2*pi/4)  exp(-j*2*pi/4)],[0.95*exp(j*2*pi/4) 0.95*exp(-j*2*pi/4)],1,1,256)
% LVxPlotZTransformMag([-1 j -j],[],1,1,256)
% LVxPlotZTransformMag([],[0.9],1,1,256)
% LVxPlotZTransformMag([(0.707*(1+j)) j (0.707*(-1+j)) -1 (0.707*(1-j)) -j -(0.707*(1+j))],[],1,1,256)
% 
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

rMax = abs(rMax);
rMin = abs(rMin);
 
 if rMin>rMax
     dtemp = rMin;
     rMin = rMax;
     rMax = dtemp;
 end
 
N = NSamps;

rGran = 0.01;
theta = 0:(2*pi/N):2*pi;

MinrMin = max(abs(PoleVec))+0.02;

if ((max(abs(PoleVec))<1.0) & (max(abs(PoleVec))>0.98))
    MinrMin = 1.0;
end

if rMin < MinrMin
    rMin = MinrMin;
end

if rMax==rMin
    rLim = 1;
else
rLim = fix((rMax - rMin)/rGran);
end

if rLim<1
    rLim=1;
end

r0 = rMin;

xMat = zeros(rLim,N+1);
yMat = zeros(rLim,N+1);

zTestVec = zeros(rLim,N+1);
zXform = ones(rLim,N+1);

for r = 1:1:rLim  
   ans = (r0 + rGran*(r-1))*exp(j*theta);
   zTestVec(r:r,1:N+1) = ans;
   
   xMat(r:r,1:N+1) = real(ans);
   yMat(r:r,1:N+1) = imag(ans);

   if isempty(PoleVec)
       zXform(r:r,1:N+1) = 1;
   else
        for poleCtr = 1:1:length(PoleVec)
        AnsForThisPole = ( 1./(1-PoleVec(poleCtr)*(1./zTestVec(r:r,1:N+1))));
        zXform(r:r,1:N+1) = zXform(r:r,1:N+1).*AnsForThisPole;
        end
   end

   if isempty(ZeroVec)
       zXform(r:r,1:N+1) = zXform(r:r,1:N+1)*1;
   else
        for zeroCtr = 1:1:length(ZeroVec)
            AnsForThisZero = 1- (ZeroVec(zeroCtr)*(1./zTestVec(r:r,1:N+1)));
            zXform(r:r,1:N+1) = zXform(r:r,1:N+1).*AnsForThisZero;
        end
   end
   
end

zXform = abs(zXform);


figure(9877)
plot3(xMat,yMat,zXform)
ylabel(['Imaginary'])
xlabel(['Real'])
zlabel(['Magnitude'])
grid on

function LVxPlotZXformMagCoeff(NumCoeffVec,DenomCoeffVec,rLim,Optionalr0,NSamps)
% function LVxPlotZXformMagCoeff(NumCoeffVec,DenomCoeffVec,rLim,Optionalr0,NSamps)
%
% LVxPlotZXformMagCoeff is intended for use with right-handed z-transforms,
% i.e., those whose ROC is the area outside of a circle having a radius just
% larger than the magnitude of the largest pole of the z-transform.
%
% NumCoeffVec is a row vector of numerator coefficients from a z-transfer fcn.
% DenomCoeffVec is a row vector of coefficients from the denominator of a
% z-transfer fcn.
% rLim is the number of circular contours in the z-plane along which to evaluate the z-transform.
% Optionalr0 is the radius of the first contour-if passed as [], the default is 1.0 for 
% FIRs (i.e., DenomCoeffVec=[1]) or the magnitude of the largest pole plus 0.02.
% If the largest pole magnitude lies betweeen 0.98 and 1.0, the radius of
% the first contour is set at 1.0 unless Optionalr0 is larger than 1.0.
% If there are no poles, pass DenomCoeffVec as [1].
% If there are no zeros, pass NumCoeffVec as [1].
% Both DenomCoeffVec and NumCoeffVec should be passed in the order of
% Coeff for z^0, Coeff for z^-1, Coeff for z^-2, etc.
% NSamps is the number of frequencies at which to evaluate the z-transform along each contour.
%
% Sample calls: 
%
% LVxPlotZXformMagCoeff([1],[1 -1.8 0.81],[1],[1],256)
% LVxPlotZXformMagCoeff([1],[1 0 0.81],[1],[1],256)
% LVxPlotZXformMagCoeff([1 1 1 1 1 1 1 1],[1],[1],[1],256)
% LVxPlotZXformMagCoeff(fir1(8,0.5),[1],[1],1,256)
% LVxPlotZXformMagCoeff(fir1(148,[0.4 0.6]),[1],[1],1,256)
% LVxPlotZXformMagCoeff([1 1 1 1 1 1 1 1],[1],[30],[0.55],256)
%
% Note:  in order to remain useful, the Polar Plot works only when rLim=1;
% You can evaluate the polar response along different radii by changing Optionalr0
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

Ang = 45;
AxFac = 1;
N = NSamps;
% optimum for average case
theta = 0:(2*pi/N):2*pi;
rGran = 0.02; % 0.02 is standard 

PoleVec = roots(DenomCoeffVec);
if isempty(PoleVec)
    PoleVec = 0;
end
MinrMin = max(abs(PoleVec))+ 0.02;

if ((max(abs(PoleVec))<1.0) & (max(abs(PoleVec))>0.98))
    MinrMin = 1.0
end

if ~(isempty(Optionalr0))
        if Optionalr0 < MinrMin
            r0 = MinrMin;
        else
            r0 = Optionalr0; 
        end
else
		if DenomCoeffVec==[1] 
   		r0 = 1.0; 
		else
		r0 = MinrMin;
        end
end

xMat = zeros(rLim,N+1);
yMat = zeros(rLim,N+1);

zTestVec = zeros(rLim,N+1);
zXform = zeros(rLim,N+1);

for r = 1:1:rLim 
 ans = (r0 + rGran*(r-1))*exp(j*theta);     
 zTestVec(r:r,1:N+1) = ans; 
 xMat(r:r,1:N+1) = real(ans);
 yMat(r:r,1:N+1) = imag(ans);  
den(1,1:N+1) = zeros(1,N+1);
nzCoef = find(abs(DenomCoeffVec)>0);

for Ctr = 1:1:length(nzCoef)
AnsForThisDCoeff = DenomCoeffVec(nzCoef(Ctr))*(zTestVec(r:r,1:N+1).^(-nzCoef(Ctr)+1));    
den(1,1:N+1) = den(1,1:N+1) + AnsForThisDCoeff;
end

num(1,1:N+1) = zeros(1,N+1);
nzCoef = find(abs(NumCoeffVec)>0);

for Ctr = 1:1:length(nzCoef)
AnsForThisNCoeff = NumCoeffVec(nzCoef(Ctr))*(zTestVec(r:r,1:N+1).^(-nzCoef(Ctr)+1));   
num(1,1:N+1) = num(1,1:N+1) + AnsForThisNCoeff;
end

zXform(r:r,1:N+1) = num(1,1:N+1)./den(1,1:N+1);
phzXform = unwrap(angle(zXform'))';
zXform(r:r,1:N+1) = abs(zXform(r:r,1:N+1));

end

figure(98762)
szZXf = size(zXform);
plot3(xMat,yMat,zXform)
grid on
ylabel(['Imaginary'])
xlabel(['Real'])
zlabel(['Magnitude'])

if rLim>1
   return
end

return
figure(98761)
grid off
zXform = zXform/(max(zXform));
xPlot = xMat.*zXform;
yPlot = yMat.*zXform;
plot(xPlot,yPlot)
xlabel(['Polar Plot of Mag(z-Transform); X-Axis = Amplitude; Contour Radius = ',num2str(r0)])
ylabel(['Amplitude'])
axis([-1.1,1.1,-1.1,1.1])




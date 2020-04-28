function MLPlotZXformMagCoeff(NumCoeffVec,DenomCoeffVec,rLim,Optionalr0,NSamps)
% function MLPlotZXformMagCoeff(NumCoeffVec,DenomCoeffVec,rLim,Optionalr0,NSamps)
%
% MLPlotZXformMagCoeff is intended for use with right-handed z-transforms,
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
% MLPlotZXformMagCoeff([1],[1 -1.8 0.81],[1],[1],256)
% MLPlotZXformMagCoeff([1],[1 0 0.81],[1],[1],256)
% MLPlotZXformMagCoeff([1 1 1 1 1 1 1 1],[1],[1],[1],256)
% MLPlotZXformMagCoeff(fir1(8,0.5),[1],[1],1,256)
% MLPlotZXformMagCoeff(fir1(148,[0.4 0.6]),[1],[1],1,256)
% MLPlotZXformMagCoeff([1 1 1 1 1 1 1 1],[1],[30],[0.55],256)
%
% Note:  in order to remain useful, the Polar Plot works only when rLim=1;
% You can evaluate the polar response along different radii by changing Optionalr0
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

Ang = 45;
AxFac = 1;
N = NSamps;
ftsz1 = 16;
ftsz2 = 14;
% optimum for average case
theta = 0:(2*pi/N):2*pi;
rGran = 0.02; % 0.02 is standard 

PoleVec = roots(DenomCoeffVec);
MinrMin = max(abs(PoleVec))+0.02;

if ((max(abs(PoleVec))<1.0) & (max(abs(PoleVec))>0.98))
    MinrMin = 1.0;
end

if ~isempty(Optionalr0)
        if Optionalr0 < MinrMin
            r0 = MinrMin;
        else
            r0 = Optionalr0; % first contour has radius Optionalr0
        end
else
		if DenomCoeffVec==[1] % no poles
   		r0 = 1.0; % it's an FIR, we'll start here since there are many poles at the origin 
   			 		 % which make the z-transform magnitude "blow-up" toward the origin
		else
		r0 = MinrMin;
        end
end

% theValOfr0 = r0

xMat(1:rLim,1:N+1) = zeros(rLim,N+1);
yMat(1:rLim,1:N+1) = zeros(rLim,N+1);

zTestVec(1:rLim,1:N+1) = zeros(rLim,N+1);
zXform(1:rLim,1:N+1) = zeros(rLim,N+1);

for r = 1:1:rLim 
   ans = (r0 + rGran*(r-1))*exp(j*theta); % z-transform is computed at multiple frequencies between 0 and 2*pi
   % lying on a constant radius contour for each value of r.     
   zTestVec(r:r,1:N+1) = ans;
   
   xMat(r:r,1:N+1) = real(ans);
   yMat(r:r,1:N+1) = imag(ans);
   
den(1,1:N+1) = zeros(1,N+1);
nzCoef = find(abs(DenomCoeffVec)>0);
for Ctr = 1:1:length(nzCoef)
AnsForThisDCoeff = DenomCoeffVec(nzCoef(Ctr))*(zTestVec(r:r,1:N+1).^(-nzCoef(Ctr)+1));    
den(1,1:N+1) = den(1,1:N+1) + AnsForThisDCoeff;
end
 
% code below works but is inefficient for large, sparse coefficient vectors 
%   for DCoefCtr = 0:1:length(DenomCoeffVec)-1 
%      AnsForThisDCoeff = DenomCoeffVec(DCoefCtr+1)*zTestVec(r:r,1:N+1).^(-DCoefCtr);
%      zXform(r:r,1:N+1) = zXform(r:r,1:N+1) + AnsForThisDCoeff;
%   end
%   zXform(r:r,1:N+1) = 1./zXform(r:r,1:N+1);

num(1,1:N+1) = zeros(1,N+1);
nzCoef = find(abs(NumCoeffVec)>0);
for Ctr = 1:1:length(nzCoef)
AnsForThisNCoeff = NumCoeffVec(nzCoef(Ctr))*(zTestVec(r:r,1:N+1).^(-nzCoef(Ctr)+1));   
num(1,1:N+1) = num(1,1:N+1) + AnsForThisNCoeff;
end

%   for NCoefCtr = 0:1:length(NumCoeffVec)-1
%      AnsForThisNCoeff = NumCoeffVec(NCoefCtr+1)*(zTestVec(r:r,1:N+1)).^(-NCoefCtr);
%      zXformZ(r:r,1:N+1) = zXformZ(r:r,1:N+1)  + AnsForThisNCoeff;
%   end
   
%   zXform(r:r,1:N+1) =  zXform(r:r,1:N+1).* zXformZ(r:r,1:N+1); % combines denom and num

zXform(r:r,1:N+1) = num(1,1:N+1)./den(1,1:N+1);

phzXform = unwrap(angle(zXform'))';
zXform(r:r,1:N+1) = abs(zXform(r:r,1:N+1));

end

a = max(max(abs(zXform)));
b = min(min(abs(zXform)));

figure(98762)
clf
set(98762,'color',[1,1,1])
hold on
szZXf = size(zXform);
plot3(xMat, yMat,zeros(szZXf(1),szZXf(2)),'b')
plot3(xMat,yMat,zXform,'b')
hold off
set(gca,'fontsize',ftsz2)
ylabel(['Imaginary Axis'],'fontsize',ftsz1)
xlabel(['Real Axis'],'fontsize',ftsz1)
zlabel(['Magnitude'],'fontsize',ftsz1)
pl = 1.1*max(abs(rLim));
axis([-pl, pl, -pl, pl, 0, 1.1*a])

grid on


if rLim>1
   return
end

figure(98761)
clf
set(98761,'color',[1,1,1])
hold on

LenzXform = length(zXform);
zXform = zXform/(max(zXform));
FreqZeroMax = 1.1*zXform(1);
FreqHalfBandMax = 1.3*max( zXform( fix(LenzXform/8):fix(LenzXform/2)-1 ) ); 
FreqNyqMax = -1.1*max(zXform( (fix(LenzXform/4)+1):fix(LenzXform/2) ) );

FreqNyqMax = (FreqNyqMax - FreqZeroMax)/2;
HorMax = max([abs(FreqNyqMax)  abs(FreqZeroMax)]);

xPlot = xMat.*zXform;
yPlot = yMat.*zXform;
plot(xPlot,yPlot,'b')
magPlot = max(sqrt(xPlot.^2 + yPlot.^2));
set(gca,'fontsize',ftsz2)
xlabel(['Polar Plot of Mag(z-Transform); X-Axis = Amplitude; Contour Radius = ',num2str(r0)],'fontsize',ftsz2)
ylabel(['Amplitude'],'fontsize',ftsz2)

theLargest = max([ max(abs(xPlot))  max(abs(yPlot)) ]);
XAxMax = 1.2*theLargest;
XAxMin = -1.2*theLargest;
YAxMin = -1.2*theLargest;
YAxMax = 1.2*theLargest;

line([XAxMin XAxMax],[0 0],'color','k','linestyle',':');
%set(hndXAxis,'color','k','linestyle',':')
line([0 0],[YAxMin YAxMax],'color','k','linestyle',':');
%set(hndYAxis,'color','k','linestyle',':')
axis([XAxMin, XAxMax, YAxMin, YAxMax])


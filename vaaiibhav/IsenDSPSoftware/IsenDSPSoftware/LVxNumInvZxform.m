function LVxNumInvZxform(NumCoefVec,DenCoefVec,M,ContourRad,nVals)
% function LVxNumInvZxform(NumCoefVec,DenCoefVec,M,ContourRad,nVals)
%
% Intended to perform the inverse z-transform for transforms whose ROC lies
% outside a circle of radius equal to the largest pole magnitude of the
% z-transform (i.e., causal sequences)
% NumCoefVec is the numerator coefficient vector or b coefficients of the
% z-transform
%
% DenCoefVec is the denominator coefficient vector or a coefficients
%
% M is the number of samples to use along the contour
%
% ContourRad is the radius of the circle used as the contour
%
% nVals is a vector of time domain indices to compute, i.e., x[nVals] is
% computed via the inverse z-transform. nVals must be a vector of
% non-negative integers corresponding to the sample indices of a causal
% sequence x[n]
%
% Sample Call:
%
% LVxNumInvZxform([1],[1 -0.9],10000,1,[0:1:20])
% LVxNumInvZxform([0 1],[1  -2  1],10000,1.05,[0:1:50])
% LVxNumInvZxform([0 1],[1 -2 1],10000,1.1,[0:1:20])
% LVxNumInvZxform([1],[1  -1.2  0.81],10000,0.91,[0:1:100])
% LVxNumInvZxform([1],[1  -1.2  0.81],10000,1,[0:1:100])
%
% LVxNumInvZxform([1,1,1,1],[1],10000,1,[0:1:20])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
ftsz1 = 16;

dw = 2*pi/M;
RadArg = [0:dw:2*pi];
ZZ = ContourRad*exp(j*RadArg);
dC = ZZ(1,2:length(ZZ)) - ZZ(1,1:length(ZZ)-1);
%ZZ = ZZ(1,1:length(ZZ)-1); % freq 0:dw:2*pi-dw

N = length(ZZ);
%compute X(z) first 
den = zeros(1,N);
nzCoeffInd = find(abs(DenCoefVec)>0);
for Ctr = 1:1:length(nzCoeffInd)
AnsForThisDCoeff = DenCoefVec(nzCoeffInd(Ctr))*(ZZ.^(-nzCoeffInd(Ctr)+1));  
den = den + AnsForThisDCoeff;
end

num = zeros(1,N);
nzCoeffInd = find(abs(NumCoefVec)>0);
for Ctr = 1:1:length(nzCoeffInd)
AnsForThisNCoeff = NumCoefVec(nzCoeffInd(Ctr))*(ZZ.^(-nzCoeffInd(Ctr)+1));   
num = num + AnsForThisNCoeff;
end

zXform = num./den;
%szzXform = size(zXform)
avgzXform = 0.5*(zXform(1,2:length(zXform))  +  zXform(1,1:length(zXform)-1));
%szavgzXform = size(avgzXform)
avgZZ = 0.5*(ZZ(1,2:length(ZZ)) + ZZ(1,1:length(ZZ)-1));

x = zeros(1,length(nVals));
% do contour integral 
for n = nVals
% x(n+1) = (1/(2*pi*j))*sum(((ZZ.^(n-1)).*zXform).*dC); % rectangular,
            % converges slowly, needs 60,000 points
x(n+1) = (1/(2*pi*j))*sum(((avgZZ.^(n-1)).*avgzXform).*dC); % averaged zXform values 
% (trapezoidal rather than rectangular)--okay with about 3000-5000 points
% for 4-5 decimal places of accuracy
end
x = x(nVals+1)

figure(120) % suitable when you know beforehand that x[n] is real
clf

xvec = nVals;
hold on
plot(xvec,real(x),'bo');
for ctr = 1:1:length(x)
    line([xvec(ctr)  xvec(ctr)],[0  real(x(ctr))])
end
ylabel('Real(x[n])')
xlabel('n')
hold off

%figure(119)

%stem3(nVals,real(x),imag(x))  
%xlabel('n')
%ylabel('Real(x[n])')
%zlabel('Imag(x[n])')
%axis([-inf inf min(real(x))  max(real(x))  min(real(x))  max(real(x))  ])
%axis([-inf inf  -inf inf  min( [-0.05 min(imag(x))]) max([0.3 max(imag(x))])     ])
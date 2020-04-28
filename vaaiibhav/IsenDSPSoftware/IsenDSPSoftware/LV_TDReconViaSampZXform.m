function LV_TDReconViaSampZXform(b,a,NumzSamps,n)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
% LV_TDReconViaSampZXform([ones(1,4)],[1],8,[-15:1:15])
% LV_TDReconViaSampZXform([ones(1,4)],[1],5,[-15:1:15])
% LV_TDReconViaSampZXform([ones(1,4)],[1],3,[-15:1:15])
kz = 0:1:NumzSamps-1;
z = exp(j*2*pi*(kz/NumzSamps));
Nm = 0; Dm = 0;
for nn = 0:-1:-length(b)+1
Nm = Nm + (b(-nn+1))*(z.^nn); end
for d = 0:-1:-length(a)+1
Dm = Dm + (a(-d+1))*(z.^d); end
zSamps = Nm./Dm;
figure(125); 
clf; N = length(zSamps);
W = exp(j*2*pi/N); k = 0:1:N-1; 
IDFS = real((W.^((n')*k))*conj(zSamps')/N);
stem(n,IDFS)


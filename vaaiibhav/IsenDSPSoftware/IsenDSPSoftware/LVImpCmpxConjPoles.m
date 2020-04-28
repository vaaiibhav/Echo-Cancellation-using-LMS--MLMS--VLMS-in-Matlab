function LVImpCmpxConjPoles(P,N)
% LVImpCmpxConjPoles(0.9,24)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    cP = conj(P); Imp = zeros(1,N);
    for n = 0:1:N; cVal = 0;
    for m = 0:1:n
     cVal = cVal + (P.^(n-m)).*(cP.^m);
    end
    if abs(imag(cVal))<10^(-15)
    cVal = real(cVal); end
    Imp(1,n+1) = cVal; end
    figure(8); subplot(211); stem(Imp)
    testImp = [1, zeros(1,N)];
    y = filter(1,[1, -P],testImp);
    y = filter(1,[1, -cP],y);
    if abs(imag(y))<10^(-15)
    y = real(y); end
    subplot(212); stem(y)
function [CCPrs] = LVCmplxConjOrd(x,tol)
% [CCPrs] = LVCmplxConjOrd(x,tol)
% x is a vector of complex conjugate poles and real poles, in radom order
% tol defines how close a pole must be to its conjugate to be detected as
% such (tol is necessary due to roundoff error, which may prevent true
% conjugates from being so detected)
% CCPrs is a vector of pairs of complex conjugates, with real poles at the
% trailing end.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Lenx = length(x);
CCPrs = x;
for ctr = 1:2:Lenx-1
   absdiff = abs(CCPrs(ctr)-conj(CCPrs(ctr+1:Lenx)));
    y = find( absdiff  < tol );
    if isempty(y) % real number
        temp = CCPrs(ctr);
        CCPrs(ctr:Lenx-1) = CCPrs(ctr+1:Lenx);
        CCPrs(Lenx) = temp;
    else
        temp = CCPrs(ctr+1);
        CCPrs(ctr+1) = CCPrs(y(1)+ctr);
        CCPrs(y(1)+ctr) = temp;
    end
end

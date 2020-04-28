function [Bc,Ac,Gain] = LVDirToCascade(b,a)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    B0 = b(1); A0 = a(1); a = a/A0; b = b/B0;
    Gain = B0/A0; lenB = length(b); lenA = length(a);
    if lenB<lenA
        b = [b,zeros(1,lenA-lenB)];
    elseif lenA<lenB
        a = [a,zeros(1,lenB-lenA)];
    end
    Ord = length(roots(b)); 
    if Ord-2*fix(Ord/2) > 0; b = [b,0]; a = [a,0]; end   
    [NumRtsB] = LVCmplxConjOrd(roots(b),10^(-10));
    [NumRtsA] = LVCmplxConjOrd(roots(a),10^(-10));
    LenNumRts = length(NumRtsB); NoPrs = fix(LenNumRts/2)
    for Ctr = 1:1:NoPrs; ind = [2*Ctr-1, 2*Ctr];
    Bc(Ctr,1:3) = real(poly(NumRtsB(ind))); 
    Ac(Ctr,1:3) = real(poly(NumRtsA(ind)));
    end

    
    


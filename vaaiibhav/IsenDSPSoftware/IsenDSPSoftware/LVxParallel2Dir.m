function [b,a] = LVxParallel2Dir(Bp,Ap,Cp)
% function [b,a] = LVxParallel2Dir(Bp,Ap,Cp)
% Receives a set of Parallel Form coefficients Bp, Ap, and Cp, and converts
% them to a set of Direct Form coefficients [b,a]
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
RRs = []; PPs = [];
szA = size(Ap);
nd = szA(1);
    
lim = nd;

for ind = 1:1:lim
    B_curr = Bp(ind,:);
    A_curr = Ap(ind,:);
    if B_curr(1,2)==0
       B_curr = B_curr(1,1);
    end  
    if A_curr(1,3)==0
       A_curr = A_curr(1,1:2);
    end  
    [r,p,k] = residuez(B_curr,A_curr);
    r = r(:)';
    p = p(:)';
    RRs = [RRs,r];
    PPs = [PPs,p];    
end

[b,a] = residuez(RRs,PPs,Cp);

ImB = imag(b);
ImB(find( abs(ImB)<10^(-14) )) = 0;
b = real(b) + j*ImB;

ImA = imag(a);
ImA(find( abs(ImA)<10^(-14) )) = 0;
a = real(a) + j*ImA;



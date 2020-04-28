function [Bp,Ap,Cp] = LVxDir2Parallel(b,a)
% Receives a set of Direct Form coefficients b and a and computes the
% Parallel Form coefficients Bp, Ap, and Cp
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
[R,P,Cp] = residuez(b,a);

P_CCOrd = LVCmplxConjOrd(P,10^(-10));
% original R and P match, so find where in P each pole of P_CCOrd is, then
% reorder R in the same manner.
[Rindices] = CCPoleLoc(P,P_CCOrd);
R_CCOrd = R(Rindices);
% proceed in pairs to collect poles and residues, and use residuez in
% reverse to obtain a and b coefficients for each second order section
lenP = length(P_CCOrd);
if rem(lenP,2)==0 % even number of poles 
    limInd = lenP-1;
    Ap = zeros(floor(lenP/2),3);
    Bp = zeros(floor(lenP/2),2);
else
    limInd = lenP-2;
    Ap = zeros(floor(lenP/2)+1,3);
    Bp = zeros(floor(lenP/2)+1,2);
end

for ind = 1:2:limInd;
p_curr = P_CCOrd(ind:ind+1);
r_curr = R_CCOrd(ind:ind+1);
[B_curr,A_curr] = residuez(r_curr,p_curr,[0]);
B_curr = B_curr(1,1:2);
%[B_curr,A_curr] = residuez(r_curr,p_curr,[])
Ap(fix(ind+1)/2,:) = real(A_curr);
Bp(fix(ind+1)/2,:) = real(B_curr);
end

if ~(rem(lenP,2)==0) % odd number of poles
   p_curr = P_CCOrd(lenP);
    r_curr = R_CCOrd(lenP);
[B_curr,A_curr] = residuez(r_curr,p_curr,[0]);
B_curr = B_curr(1,1:1);
%[B_curr,A_curr] = residuez(r_curr,p_curr,[]);
 Ap((lenP+1)/2,:) = [real(A_curr),0];
Bp((lenP+1)/2,:) = [real(B_curr),0];   
end

function [Rindices] = CCPoleLoc(P,P_CCOrd)
Rindices = [];
for PCC_Ctr = 1:1:length(P_CCOrd)
    for PCtr = 1:1:length(P)
        if abs(P(PCtr)-P_CCOrd(PCC_Ctr)) < 10^(-10)
            Rindices = [Rindices,PCtr];
        end
    end
end



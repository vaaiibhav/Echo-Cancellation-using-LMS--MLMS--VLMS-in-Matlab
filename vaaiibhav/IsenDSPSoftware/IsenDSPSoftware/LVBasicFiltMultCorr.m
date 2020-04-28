function Imp = LVBasicFiltMultCorr(N,LoK,HiK)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
% Imp = LVBasicFiltMultCorr(30,0,10)
Imp = 0; 
for k = LoK:1:HiK
if k==0|k==N/2
A = 1; else
A = 2; end
Imp = Imp + A*((-1)^k)*cos( 2*pi*k*( 0:1:N-1 )/ N );
end

lenImp = length(Imp)
LVFreqResp(Imp, 500)
function LVDesignDigFiltViaImpInv(Wp,Ws,Rp,As,FiltType,Fs)
% LVDesignDigFiltViaImpInv(0.4*pi,0.5*pi,1,40,1,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
OmP = Wp*Fs; OmS = Ws*Fs;
if FiltType==1
[Z,P,K] = LVDesignButterworth(OmP,OmS,Rp,As);
elseif FiltType==2
 %   [Z,P,K] = LVDesignCheby1BPF(4,5,8,9,Rp,As)   
[Z,P,K] = LVDesignCheby1Filter(Rp,As,OmP,OmS);
elseif FiltType==3
[Z,P,K] = LVDesignCheb2(Rp,As,OmP,OmS);
else
[Z,P,K] = LVDesignEllip(Rp,As,OmP,OmS);
end; b = real(K*poly(Z)), a = real(poly(P)),
[BZ,AZ] = impinvar(b,a,Fs);
LVsFRzFrLog(b,a,OmP,BZ,AZ,97)

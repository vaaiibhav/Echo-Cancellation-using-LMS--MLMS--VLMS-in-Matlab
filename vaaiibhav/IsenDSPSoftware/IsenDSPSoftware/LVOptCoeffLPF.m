function LVOptCoeffLPF(L,wp,ws,THi,TLo,Dec)  
% LVOptCoeffLPF(40,0.5,0.6,0.5,0.35,0.01)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Dec = -abs(Dec); noComp = ceil((THi-TLo)/abs(Dec)); ctr = 0;
TSB = zeros(noComp,3); LenFFT = 4096; if rem(L,2)==0; limK = L/2-1; else;
limK = (L-1)/2; end; n = 0:1:L-1; M=(L-1)/2;
for T = THi:Dec:TLo; Ak = [ones(1,11), T, zeros(1,8)]; 
LA=length(Ak); WF=(cos(((n-M)')*[0:1:LA-1]*2*pi/L))*([Ak(1),2*Ak(2:LA)]');
WF = WF'; Imp = WF/L; figure(3); 
subplot(211); stem([0:1:length(Imp)-1],Imp); xlabel('n'); ylabel('Amplitude'); 
fr = abs(fft(Imp,LenFFT)); fr = fr(1,1:LenFFT/2+1); fr = fr/(max(abs(fr))); 
Lfr = length(fr); SB = fr(1,round(ws*(Lfr-1)):(Lfr-1)); 
PB = fr(1,1:round((wp/pi)*(Lfr-1))); subplot(212); 
SBAt = -20*log10(max(SB) +eps); PBR = -20*log10(min(PB)); 
ctr = ctr + 1; TSB(ctr,1)=T; TSB(ctr,2)=SBAt; TSB(ctr,3)=PBR;
plot([0:1:Lfr-1]/(Lfr-1),20*log10(fr+eps)); strSBA = num2str(SBAt);
xlabel(['Norm.Freq.(T = ',num2str(T),' and As = ',strSBA,')'])
ylabel('Magnitude'); axis([0 1 -70 10]); pause(0.01); 
end; sc = TSB(:,2); bestSB = max(sc), bestSBind = find(sc==bestSB); 
bestT = TSB(bestSBind,1), finalPBRipple = TSB(bestSBind,3)
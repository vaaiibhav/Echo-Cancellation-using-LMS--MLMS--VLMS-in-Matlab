function LVButterFR(N,OmegaC)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
% LVButterFR(3,1)
k = 0:1:2*N-1; 
P = OmegaC*exp(j*(pi/(2*N))*(2*k+N+1));
NetP = P(find(real(P)<0))
Freq = (0:0.01:8*OmegaC); s = j*Freq;
Resp = OmegaC^N; for Ctr = 1:1:length(NetP)
Resp = Resp.*(1./(s - NetP(Ctr)));
end; figure(97); subplot(311); plot(Freq,20*log10(abs(Resp)));
xlabel('(a) Freq, Rad/s'); ylabel('Mag, dB')
subplot(312); plot(Freq,abs(Resp)); xlabel('(b) Freq, Rad/s'); 
ylabel('Mag'); subplot(313); plot(Freq,unwrap(angle(Resp)))
xlabel('(c) Freq, Rad/s'); ylabel('Radians')
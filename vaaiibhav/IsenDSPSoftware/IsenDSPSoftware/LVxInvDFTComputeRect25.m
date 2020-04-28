function LVxInvDFTComputeRect25
% LVxInvDFTComputeRect25
% Demonstrates the computation of the Inverse 
% Discrete Fourier Transform (DFT) using a direct 
% implementation. After calling the function, press any key for the next
% computation.
% The call:
% LVxInvDFTComputeRect25
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
N = 32;
t = 0:1/32:1-1/32;  % = (0:1:N-1)/N;

Waveform = [ones(1,8), -ones(1,24)];
% This is a rectangular wave of 25% duty cycle, amplitudes
% balanced around zero.
x = 0:1:31;
% First, compute the DFT of the test waveform to use 
% in synthesis

CumInvDFT = zeros(1,length(t));
figure(68)

subplot(321)
stem(x,Waveform);
ylabel(['Amplitude'])
xlabel(['(a)  The Original Signal, X=Sample'])
axis([0, 32, -1.5, 1.5])

A = 1/32;

for harmno = 0:1:31 % this part obtains the DFT bin coefficients so they can be used below
testimag = -sin(2*pi*t*(harmno));
% roundoff error at harmno = 16 for sine component can cause very tiny
% non-zero imag component at bin 16
    if abs(testimag<(10^-12)) 
        testimag=0;
    end
testreal = cos(2*pi*t*(harmno));
DFTrealbin(1,harmno+1) = A*sum(Waveform.*testreal);
DFTimagbin(1,harmno+1) = A*sum(Waveform.*testimag);
end
% DFTCoeff = fft(Waveform); % the simple way using Matlab, replaces above
% loop
DFTCoeff = DFTrealbin + i*DFTimagbin;  % Bins 0 to 31
% Synthesis portion--the IDFT itself

for harmno = 0:1:16 % includes Bins 0, N/2 (16) and all positive; 
                    % negative freqs are Bins 17-31 which are equivalent to Bins -15:-1:-1
SynthExp = exp(j*2*pi*t*(harmno)); % = cos(2*pi*t*(harmno)) + i*sin(2*pi*t*(harmno))
IDFTPosK = DFTCoeff(1,harmno+1)*SynthExp; 
IDFTNegK = conj(DFTCoeff(1,harmno+1))*conj(SynthExp); % takes advantage of conjugate symmetry, 
% based on signal being real

if harmno==0|harmno==16
CumInvDFT = CumInvDFT + IDFTPosK; 
else
CumInvDFT = CumInvDFT + IDFTPosK + IDFTNegK;     
end

subplot(321)
stem(x,Waveform);
ylabel(['Amplitude'])
xlabel(['(a)  The Original Signal, X=Sample'])
axis([0, 32, -1.5, 1.5])

subplot(323)
theRealOut = real(IDFTPosK);
stem(x,theRealOut);
ylabel(['Real'])
xlabel(['(c) Real(IFFT), Bin ',num2str(harmno),'; X=Sample'])
axis([0, 31, -0.6, 0.6])

subplot(324)
theImagOut = imag(IDFTPosK);
stem(x,theImagOut);
ylabel(['Imag'])
xlabel(['(d) Imag(IFFT), Bin ',num2str(harmno),'; X=Sample'])
axis([0, 31, -0.6, 0.6])

%--------------------
subplot(325)
if ~(harmno==0|harmno==16)
theRealOut = real(IDFTNegK);
stem(x,theRealOut);
ylabel(['Real'])
xlabel(['(e) Real(IFFT), Bin -',num2str(harmno),'; X=Sample'])
else
stem(x,zeros(1,length(IDFTNegK)),'bo');
ylabel(['Real'])
xlabel(['(e) No Neg Component for Bin ',num2str(harmno)])
end
axis([0, 31, -0.6, 0.6])

subplot(326)
if ~(harmno==0|harmno==16)
theImagOut = imag(IDFTNegK);
stem(x,theImagOut);
ylabel(['Imag'])
xlabel(['(f) Imag(IFFT), Bin -',num2str(harmno),'; X=Sample'])
else
stem(x,zeros(1,length(IDFTNegK)));
ylabel(['Imag'])
xlabel(['(e) No Neg Component for Bin ',num2str(harmno)])
end
axis([0, 31, -0.6, 0.6])

subplot(322)
if ~(harmno==0|harmno==16)
theRealCum = CumInvDFT;
stem(t*32,theRealCum);
ylabel(['Amplitude'])
xlabel(['(b) Sum(IFFT) to Bins +/- ',num2str(harmno),'; X=Sample'])
else
theRealCum = CumInvDFT;
stem(t*32,theRealCum ,'bo');
ylabel(['Amplitude'])
xlabel(['(b) Sum(IFFT) to Bin ',num2str(harmno),'; X=Sample'])    
end
axis([0, 32, -1.5, 1.5])
if harmno==16
   return
end
  
pause

end

 

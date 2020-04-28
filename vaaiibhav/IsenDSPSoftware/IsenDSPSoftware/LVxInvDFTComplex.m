function LVxInvDFTComplex(theComplexSignal)
% function LVxInvDFTComplex(theComplexSignal)
% Demonstrates the computation of the Inverse 
% Discrete Fourier Transform (DFT) for a complex signal input by the user.
% Test calls:
% LVxInvDFTComplex(([(1+j),(2-j),-1,(3-2*j)]))
% LVxInvDFTComplex([exp(j*2*pi*3.5*(0:1:31)/32)])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

Waveform = theComplexSignal;

N = length(theComplexSignal);
n = 0:1:N-1;

if rem(N,2)==0 % N even
    ReconLim = N/2;
else
    ReconLim = (N-1)/2;
end

% First, compute the DFT of the test waveform to use 
% in synthesis

CumInvDFT = 0;
figure(68)

xvec = 0:1:N-1;

subplot(421)
stem(xvec,real(Waveform),'bo');
ylabel(['Real'])
xlabel(['(a)  The Original Signal, X=Sample'])
plotlim = 1.2*max(abs(real(Waveform)));
axis([0, N, -plotlim, plotlim])

subplot(422)
stem(xvec,imag(Waveform),'bo');
ylabel(['Imag'])
xlabel(['(b)  The Original Signal, X=Sample'])
axis([0, N, -inf, inf])

A = 1/N;

for harmno = 0:1:N-1 % this part obtains the DFT bin coefficients so they can be used below
testimag = -sin(2*pi*n/N*(harmno));

testreal = cos(2*pi*n/N*(harmno));
DFTrealbin(1,harmno+1) = A*sum(Waveform.*testreal);
DFTimagbin(1,harmno+1) = A*sum(Waveform.*testimag);
end

% Synthesis portion--the IDFT itself
theDFT = DFTrealbin + j*DFTimagbin;

for harmno = 0:1:ReconLim
    PosK = harmno;
    NegK = N - harmno;
    
IDFTPosK = theDFT(1,PosK+1)*exp(j*2*pi*n*PosK/N); 

if ~(harmno==0|harmno==N/2)
IDFTNegK = theDFT(1,NegK+1)*exp(j*2*pi*n*NegK/N);  
end

if harmno==0|harmno==(N/2)
CumInvDFT = CumInvDFT + IDFTPosK; 
else
CumInvDFT = CumInvDFT + IDFTPosK + IDFTNegK;     
end

subplot(425)
theRealOut = real(IDFTPosK);
stem(xvec,theRealOut,'bo');
ylabel(['Real(IFFT)'])
xlabel(['(e) Real(IFFT), Bin ',num2str(harmno),'; X=Sample'])
axis([0, N, -inf, inf])

subplot(426)
theImagOut = imag(IDFTPosK);
stem(xvec,theImagOut,'bo');
ylabel(['Imag(IFFT)'])
xlabel(['(f) Imag(IFFT), Bin ',num2str(harmno),'; X=Sample'])
axis([0, N, -inf, inf])

%--------------------
subplot(427)
if ~(harmno==0|harmno==(N/2))
theRealOut = real(IDFTNegK);
stem(xvec,theRealOut,'bo');
ylabel(['Real(IFFT)'])
xlabel(['(g) Real(IFFT), Bin -',num2str(harmno),'; X=Sample'])
else
   stem(xvec,zeros(1,length(theRealOut)),'bo');
   xlabel(['(g) No Neg Component for Bin ',num2str(harmno)])
   ylabel([''])
end
axis([0 N -inf inf])

subplot(428)
if ~(harmno==0|harmno==(N/2))
theImagOut = imag(IDFTNegK);
stem(xvec,theImagOut,'bo');
ylabel(['Imag(IFFT)'])
xlabel(['(h) Imag(IFFT), Bin -',num2str(harmno),'; X=Sample'])
else
  stem(xvec,zeros(1,length(theImagOut)),'bo');
  xlabel(['(h) No Neg Component for Bin ',num2str(harmno)])
  ylabel([''])
end
axis([0 N -inf inf])

subplot(423)
if ~(harmno==0|harmno==(N/2))
theRealCum = real(CumInvDFT);
stem(xvec,theRealCum ,'bo');
ylabel(['Real'])
xlabel(['(c) Sum(IFFT) to Bins +/- ',num2str(harmno),'; X=Sample'])
else
theRealCum = real(CumInvDFT);
stem(xvec,theRealCum ,'bo');
ylabel(['Real'])
xlabel(['(c) Sum(IFFT) to Bin ',num2str(harmno),'; X=Sample'])    
end
axis([0, N, -plotlim, plotlim])

subplot(424)
if ~(harmno==0|harmno==(N/2))
theImagCum = imag(CumInvDFT);
stem(xvec,theImagCum ,'bo');
ylabel(['Imag'])
xlabel(['(d) Sum(IFFT) to Bins +/- ',num2str(harmno),'; X=Sample'])
else
theImagCum = imag(CumInvDFT);
stem(xvec,theImagCum ,'bo');
ylabel(['Imag'])
xlabel(['(d) Sum(IFFT) to Bin ',num2str(harmno),'; X=Sample'])    
end
axis([0 N -inf inf])


if harmno==N/2
   return
end
  
pause(1)

end

 

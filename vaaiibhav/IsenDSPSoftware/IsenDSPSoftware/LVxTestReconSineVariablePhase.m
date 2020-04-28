function LVxTestReconSineVariablePhase(k1,N,PhaseDegrees)
% function LVxTestReconSineVariablePhase(k1,N,PhaseDegrees)
% Performs correlations at the zeroth lag between
% test cosine and sine waves of frequency k1, having N samples, and a
% sinusoid having a phase of PhaseDegrees, and 
% then reconstructs the original sinusoid of phase equal 
% to PhaseDegrees by using the CZL values and the test sine and cosine.
%Test calls:
% LVxTestReconSineVariablePhase(1,32,45)
% LVxTestReconSineVariablePhase(0,32,90)
% LVxTestReconSineVariablePhase(16,32,90)
% LVxTestReconSineVariablePhase(1,32,45)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

k1 = fix(k1)

PhaseRadians = (PhaseDegrees/360)*(2*pi);
t = 0:1:N-1;

TestCosCorr = cos(2*pi*t*(k1/N));  
TestSinCorr = sin(2*pi*t*k1/N);

TestSine = sin(2*pi*t*k1/N + PhaseRadians);

CorrCos = sum(TestCosCorr.*TestSine);

CorrSin = sum(TestSinCorr.*TestSine);

if k1==0|k1==N/2
ReconSine = (1/N)*(TestCosCorr*CorrCos + TestSinCorr*CorrSin);    
else
ReconSine = (2/N)*(TestCosCorr*CorrCos + TestSinCorr*CorrSin);
end

figure(14811)
clf

subplot(221)
lny1Cos = length(TestCosCorr);
stem(0:1:lny1Cos-1,TestCosCorr,'bo');
axis([0 N -1.2 1.2])
xlabel(['(a)  Sample, Cosine Correlator'])
ylabel(['Amp'])

subplot(223)
stem(0:1:lny1Cos-1,TestSinCorr,'bo');
axis([0 N -1.2 1.2])
xlabel(['(c)  Sample, Sine Correlator'])
ylabel(['Amp'])

subplot(222)
thingtoplot = TestSine;
stem(0:1:N-1,TestSine,'bo');
grid on
ylabel(['Amp'])
xlabel(['(b)  Sample, Test Sig'])
axis([0,N,min([1.2*min(thingtoplot),-1.2]),max([1.2*max(thingtoplot),1.2])])

subplot(224)
otherthingtoplot = ReconSine;
stem(0:1:N-1,ReconSine,'bo');
min1 = min([1.2*min(otherthingtoplot),-1]);
max1 = max([1.2*max(otherthingtoplot),1]);   
axis([0,N,min1,max1])
grid on
ylabel(['Amp'])
xlabel(['(d)  Sample, Recon Sig'])
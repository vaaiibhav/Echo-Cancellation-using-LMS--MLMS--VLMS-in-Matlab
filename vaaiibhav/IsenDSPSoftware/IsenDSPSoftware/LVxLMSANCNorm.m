function  LVxLMSANCNorm(PlantCoeffVec,k,Mu,freq,DVMult)
% function  LVxLMSANCNorm(PlantCoeffVec,k,Mu,freq,DVMult)
% This script implements a 10-tap LMS adaptive filter in an
% Active Noise Cancellation (ANC) application, such as cancell-
% ing noise in a duct.
% PlantCoeffVec is a row vector of 10 coefficients for the FIR that
% simulates the Plant (i.e., the Duct impulse response).
% k specifies the standard deviation of random noise to be mixed with two 
% sinusoids having frequencies freq and 3*freq and respective amplitudes of
% 1.0 and 0.5;
% Mu specifies the tap weight update damping coefficient;
% An NLMS algorithm is used; its effectiveness can be tested by giving the
% test signal various stepped-amplitude profiles with the input argument
% DVMult, which is a vector of amplitudes to impose on the test signal as a
% succession of equally spaced amplitude segments over a test signal of 15 
% times the filter length.
% A typical call:
% LVxLMSANCNorm([0,0,1,0,-0.5,0.6,0,0,-1.2,0],2,2,27,[1,2,5,8])
% A call which shows slow or poor convergence due to lack of noise 
% with the sinusoidal test signal is
% LVxLMSANCNorm([0,0,1,0,-0.5,0.6,0,0,-1.2,0],0,2,3,[1,2,5,8])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

FiltLen= 10;
SR = 150;
LenHalfTD = 75;

AmpProfile = ones(1,SR);
ld = length(DVMult);

DivPts = [1,fix(SR*[1:1:(ld-1)]/ld)];
DivPts = [DivPts,(SR+1)];
ld = ld+1;
for Ctr = 1:1:ld-1
AmpProfile(1,DivPts(Ctr):DivPts(Ctr+1)-1) = DVMult(Ctr);
end

t = 0:1/(SR-1):1;
DataVec = sin(2*pi*t*freq)+ 0.5*cos(2*pi*t*freq*3)+ k*randn(1,SR);

DataVec = AmpProfile.*DataVec;

TapWt = zeros(length(DataVec),FiltLen);
FiltOut = zeros(1,length(DataVec));
Err = zeros(1,length(DataVec));
n = 1:1:FiltLen;

IterationIndices = FiltLen + 1:1:length(DataVec);
NoIts = length(IterationIndices);

for CurPtr = IterationIndices

FiltOut(1,CurPtr) = sum((TapWt(CurPtr,n)).*DataVec(CurPtr-n));
Err(1,CurPtr)= sum(PlantCoeffVec(n).*DataVec(1,CurPtr-n)) - FiltOut(1,CurPtr);
LIMAdj = sum(DataVec(1,CurPtr-n).^2)+ eps;
TapWt(CurPtr+1,n)= TapWt(CurPtr,n) + ...
    (Mu*Err(1,CurPtr).*(DataVec(1,CurPtr-n)))/(2*LIMAdj); 
end

figure(18)
clf
hold on
for colctr = 1:1:FiltLen
    plot(TapWt(:,colctr),'k')
end
grid on
xlabel(['Iteration'])
ylabel(['Amplitude'])
plotlim = 1.2*max(max(abs(TapWt)));
axis([0,SR,-plotlim,plotlim])


figure(19)
subplot(211)
plot(DataVec)
plotlim1 = 1.2*max(abs(DataVec));
xlabel(['(a) Sample, Input Data '])
ylabel('Amplitude')
axis([0 length(t) -plotlim1 plotlim1])

subplot(212)
plot(Err)
plotlim2 = 1.2*max(abs(Err));
xlabel(['(b) Sample/Iteration of Output/Error'])
ylabel('Amplitude')
axis([0 length(t) -plotlim2 plotlim2])








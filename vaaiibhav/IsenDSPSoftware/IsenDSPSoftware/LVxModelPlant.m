function LVxModelPlant(A, B,LenLMS, NoPrZs, NoPrPs, Mu, tSig, NAmp, NoIts)
%
% function LVxModelPlant(A, LenLMS, NoPrZs, NoPrPs, Mu, tSig, NAmp, NoIts)
% A and B are a set of IIR and FIR coefficients to generate the Plant
% response
% LenLMS is the number of LMS adaptive FIR coefficients to use to model the
% IIR; it should generally be long enough to model the impulse response until it
% has mostly decayed away
% NoPrZs is the number of Prony zeros to use to model the LMS-derived
% impulse response (the converged LMS adaptive filter coefficients, which form a 
% truncated version of the Plant impulse response);
% NoPrPs is the number of Prony poles to use to model the LMS-derived impulse response;
% Mu is the standard LMS update term weight;
% tSig if passed as 0 = white noise; 1 = DC; 2 = Nyquist; 3 = Half-band;
%   4 = mixed cosines (DC, Half-band, and Nyquist)
% NAmp = amplitude of noise to add to selected tSig
% NoIts is the number of iterations to do in the adaptive process. It should
% be at least 10 times LenLMS.
%
% Example calls:
%
% LVxModelPlant([1,-0.9],[1],100, 2,2, 0.5,0,0,1000)
% LVxModelPlant([1, 0,0.81],[1],100, 3,3,0.5,1, 0,1000)
% LVxModelPlant([1,-1.3,0.81],[1],100, 3,3,0.5,1, 0,1000)
% LVxModelPlant([1,-1,0.64],[1],100, 3,3, 0.5, 1, 0,1000)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
LentSig = NoIts;
if NoIts<LentSig
    NoIts=LentSig+1;
end
if tSig==0
    x = randn(1,LentSig);
    SigType = 'Random Noise';
elseif tSig==1
 x = ones(1,LentSig);  %  DC
     SigType = 'Unit Step (DC)';
elseif tSig==2 
 x = zeros(1,LentSig);
 nn = 1:1:LentSig/2;
 x(1,2*nn ) = -ones;
 x(1,2*nn -1) = ones; % Nyquist limit frequency (fs/2)
 SigType = 'Nyquist Limit (fs/2)';
elseif tSig==3
 x = zeros(1,LentSig);
 nn = 1:1:LentSig;
 x(1,4*nn ) = ones;
 x(1,4*nn -2) = -ones; % half-band frequency (fs/4)
 x = x +  eps;
 SigType = 'Half-Band Freq (fs/4)';
elseif tSig==4 % 3 cosines      
x = cos(2*pi*0*(1:1:LentSig)/LentSig) + cos(2*pi*256*(1:1:LentSig)/LentSig)+ cos(2*pi*512*(1:1:LentSig)/LentSig);
SigType = ['Mixed Cosines'];
else
 x = randn(1,LentSig);
 SigType = 'Random Noise';
end

x = x/max(abs(x));
x = x + NAmp*randn(1,length(x));
P = filter(B,A,x); % the Plant output for the input signal x

Bmat = zeros(LenLMS,NoIts);
start = LenLMS + 1;
y = zeros(1,NoIts);
rangeZ = 1:1:LenLMS; 

for Ctr = start:1:NoIts-1
  y(1,Ctr) = sum(Bmat(:,Ctr).*x(Ctr - rangeZ + 1)');
  FiltPower = sum(  x(Ctr - rangeZ + 1).^2 ) + 0.001;
  Err(Ctr) = P(Ctr) - y(1,Ctr);     
  Bmat(1:LenLMS,Ctr + 1) = Bmat(1:LenLMS,Ctr) + 2*Mu*(1/FiltPower)*Err(Ctr).*(x(Ctr-(1:1:LenLMS) + 1)');
end

figure(888)
clf

theImp = filter(B,A,[1  zeros(1,LenLMS)]);
bCoeffEst = Bmat(:,NoIts);

subplot(313)
hold on
stem(theImp,'b*');
stem(bCoeffEst,'bo');
hold off

xlabel('Sample, Actual ImpResp (Stars) of Plant & Estimated ImpResp (Circles) Via LMS')
ylabel('Amplitude')
axis([0   length(bCoeffEst)  (min([min(bCoeffEst) min(theImp)])-0.2*abs(min([min(bCoeffEst) min(theImp)])))   (max([max(bCoeffEst) max(theImp)])+0.2*abs(max(([max(bCoeffEst) max(theImp)])))) ])
plotXlim = 100;
subplot(311)
stem(x(1,1:plotXlim));
xlabel(['Sample, Test Signal (Partial)'])
ylabel('Amplitude')
axis([0  plotXlim  min([0,1.1*min(x)])  1.1*max(x)])

subplot(312)
stem(P(1,1:plotXlim));
xlabel(['Sample, Plant Output (Partial)'])
ylabel('Amplitude')
axis([0  plotXlim  min([0,1.1*min(P)])  1.1*max(P)])

OrdPrZs = NoPrZs-1;
OrdPrPs = NoPrPs-1;

[finalB, finalA] = prony(bCoeffEst, OrdPrZs, OrdPrPs)

strfinalB = '';
for ctr = 1:1:length(finalB)
strfinalB = [strfinalB,'  ',num2str(finalB(ctr),2)];
end

strfinalA = '';
for ctr = 1:1:length(finalA)
strfinalA = [strfinalA,'  ',num2str(finalA(ctr),2)];
end

strResults = [' Prony A = [',strfinalA,' ]; Prony B = [',strfinalB,' ]']

return
xdecon = filter(finalA,finalB,[1,zeros(1,80)]);

figure(876)

stem(conv(xdecon,P))





function LV_FIRPSurfTrueMSE(PC1,PC2,Mu,c1Strt,c2Strt,NoIts,tstDType, CosFrq)
% function LV_FIRPSurfTrueMSE(PC1,PC2,Mu,c1Strt,c2Strt,NoIts,tstDType,CosFrq)
%
% This functions plots the performance surface of a 2-tap adaptive FIR in 2-dimensions
% along with the coefficient values at each iteration.
%
% PC1 and PC2 are the two Plant coefficients to be determined;
% Mu is the usual LMS update term weight;
% c1Strt and c2Strt are the initial coefficient estimates;
% NoIts is the number of iterations to perform;
% tstDType is passed as 0 for white noise, 1 for a unit step (DC), 2 for
% the Nyquist Limit frequency, 3 for the half-band frequency, and 4 for a
% pure sine wave.  
% CosFrq is the frequency for the test sine wave (pass as [] if a sine wave
% is not being used as a test signal).  
%
% Typical calls:
%
% LV_FIRPSurfTrueMSE(3,2,0.25,0,0,6,1,[])
% LV_FIRPSurfTrueMSE(1,-0.5,0.35,0,-0.2,12,0,[])
% LV_FIRPSurfTrueMSE(1,-0.5,0.35,4,-5,12,3,[])  
% LV_FIRPSurfTrueMSE(-1,3,0.25,5,-4,12,0,[])
% LV_FIRPSurfTrueMSE(1,-0.5,0.35,0,-0.2,12,0,[]) 
% LV_FIRPSurfTrueMSE(1,-0.5,0.35,0,-0.2,12,4,[125]) 
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

if isempty(CosFrq)
CosFrq = 5;
end

if tstDType==0
    DataVec = randn(1,1024);
    SigType = 'Random Noise';
elseif tstDType==1
 DataVec = ones(1,1024);  %  DC
     SigType = 'Unit Step (DC)';
elseif tstDType==2 
 DataVec = zeros(1,1024);
 nn = 1:1:512;
 DataVec(1,2*nn ) = -1;
 DataVec(1,2*nn -1) = 1; % Nyquist limit frequency (fs/2)
 SigType = 'Nyquist Limit (fs/2)';
elseif tstDType==3
 DataVec = zeros(1,1024);
 nn = 1:1:256;
 DataVec(1,4*nn ) = 1;
 DataVec(1,4*nn -2) = -1; % half-band frequency (fs/4)
 DataVec = DataVec +  eps;
 SigType = 'Half-Band Freq (fs/4)';
elseif tstDType==4 % cosine      
DataVec = cos(2*pi*CosFrq*(1:1:1024)/1024) + eps;
SigType = ['Cosine, Freq = ',num2str(CosFrq)];
else
 DataVec = randn(1,1024);
 SigType = 'Random Noise';
end

figure(8654)
subplot(2,1,1)
stem(DataVec(1,1:NoIts))
ylabel(['Amplitude'])
xlabel(['Test Data Vector x'])
axis([0 NoIts -1.2*max(abs(DataVec))  1.2*max(abs(DataVec)) ])

rat = DataVec(1,2:NoIts)./DataVec(1,1:NoIts-1);

subplot(2,1,2)
stem(rat)
ylabel(['x(n+1)/x(n)'])
xlabel(['Sample Index n for x (Test Data Vector)'])
axis([0 NoIts -1.2*max(abs(rat))  1.2*max(abs(rat)) ])

c1Lim = 2*max(  [abs(c1Strt)  abs(PC1)]  );
c2Lim = 2*max(  [abs(c2Strt)  abs(PC2)]  );

%====================================================================
Testc1 = PC1-c1Lim:c1Lim/20:PC1+c1Lim;  
Testc2 = PC2-c2Lim:c2Lim/20:PC2+c2Lim;  

theMSEIdeal = zeros( length(Testc1),length(Testc2) );

for testc1Ctr = 1:1:length(Testc1)  
theMSEIdeal(testc1Ctr,1:length(Testc2)) = (PC1-Testc1(testc1Ctr))^2 + (PC2 -Testc2).^2;
end

c1 = [PC1-c1Lim,  PC1+c1Lim]; 
%==========================================================================
c1Est = zeros(1,NoIts);
c2Est = zeros(1,NoIts);

c1Est(1,1) = c1Strt;
c2Est(1,1) = c2Strt;
%==========================================================================
figure(1869);

for MCtr = 1:1:NoIts;

Err1(MCtr) = (PC1 - c1Est(MCtr))*DataVec(MCtr+1);
Err2(MCtr) = (PC2 - c2Est(MCtr))*DataVec(MCtr);

theMSE(MCtr) = Err1(MCtr)^2 + Err2(MCtr)^2;
PowerInFilt = DataVec(MCtr+1)^2 + DataVec(MCtr)^2 + 1;
c1Est(1,MCtr+1) = c1Est(1,MCtr) + 2*Mu*Err1(MCtr)*DataVec(MCtr+1)/PowerInFilt;
c2Est(1,MCtr+1) = c2Est(1,MCtr) + 2*Mu*Err2(MCtr)*DataVec(MCtr)/PowerInFilt;

matTestc1 = ((Testc1')*ones(1,length(Testc2)))';
matTestc2 = (Testc2')*ones(1,length(Testc1));

contour(matTestc1,matTestc2,theMSEIdeal)

xlabel(['Coefficient for Tap 1;  Test Signal Type is ',SigType,'; Iteration = ',num2str(MCtr)])
ylabel(['Coefficient for Tap 2'])
axis([PC1-c1Lim  PC1+c1Lim  PC2-c2Lim  PC2+c2Lim])

hold on

plot(c1Est(1,1:MCtr),c2Est(1,1:MCtr),'.');
plot(c1Est(1,1:MCtr),c2Est(1,1:MCtr)); 
plot(PC1,PC2,'bo');

%  replot  locus to include updated point

pause(0.5)

plot(c1Est(1,1:MCtr+1),c2Est(1,1:MCtr+1),'.');
plot(c1Est(1,1:MCtr+1),c2Est(1,1:MCtr+1)); 
xlabel(['Coefficient for Tap 1;  Test Signal Type is ',SigType,'; Iteration = ',num2str(MCtr)])
ylabel(['Coefficient for Tap 2'])
axis([PC1-c1Lim  PC1+c1Lim  PC2-c2Lim  PC2+c2Lim]) 
%=======================================================================================
hold off
end






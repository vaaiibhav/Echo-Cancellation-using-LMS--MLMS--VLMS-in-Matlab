function LV_FIRLMS6Panel(PC1,PC2,Mu,c1St,c2St,NoIts,tstDType, CosFrq, NoiseAmp)
% function LV_FIRLMS6Panel(PC1,PC2,Mu,c1St,c2St,NoIts,tstDType,CosFrq,NoiseAmp)
%
% This functions plots the performance surface of a 2-tap adaptive FIR in 2-dimensions
% along with the coefficient values at each iteration.
% PC1 and PC2 are the two Plant coefficients to be determined;
% Mu is the usual LMS update term weight, c1St and c2St are the
% initial coefficient estimates;
% NoIts is the number of iterations to perform;
% tstDType is passed as 0 for white noise, 1 for a unit step (DC), 2 for
% the Nyquist Limit frequency, 3 for the half-band frequency, and 4 for a pure
% sine wave.  
% SineFrq is the frequency for the test cosine wave (pass as []
% if a cosine wave is not being used as a test signal). 
% Typical calls are:
%
% LV_FIRLMS6Panel(3,2,0.25,0,0,6,1,[],0)
% LV_FIRLMS6Panel(1,-0.5,0.35,0,-0.2,12,0,[],0)
% LV_FIRLMS6Panel(1,-0.5,0.35,4,-5,6,3,[],0)   
% LV_FIRLMS6Panel(-1,3,0.25,5,-4,12,4,25,0)
% LV_FIRLMS6Panel(-1,3,0.5,5,-4,6,4,[155],0)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if isempty(CosFrq)
   CosFrq = 25;
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
 DataVec(1,4*nn ) = -1;
 DataVec(1,4*nn -2) = 1; % half-band frequency (fs/4)
 DataVec = DataVec +  0.01*randn(1,1024);
 SigType = 'Half-Band Freq (fs/4)';
elseif tstDType==4 % cosine      
DataVec = cos(2*pi*CosFrq*([1:1:1024])/1024); % +  0.05*randn(1,1024);
SigType = 'Cosine';
else
 DataVec = randn(1,1024);
 SigType = 'Random Noise';
end

DataVec = DataVec + NoiseAmp*randn(1,length(DataVec));

c1Lim = 1.25*max(  [abs(c1St)  abs(PC1)]  );
c2Lim = 1.25*max(  [abs(c2St)  abs(PC2)]  );

%====================================================================
Testc1 = PC1-c1Lim:c1Lim/20:PC1+c1Lim;  
Testc2 = PC2-c2Lim:c2Lim/20:PC2+c2Lim;  

theMSEIdeal = zeros(length(Testc1),length(Testc2));

c1 = [PC1-c1Lim  PC1+c1Lim]; 
%=================================================================================
c1EstLMS = zeros(1,NoIts);
c2EstLMS = zeros(1,NoIts);

c1EstLMS(1,1) = c1St;
c2EstLMS(1,1) = c2St;

ErrLMS = zeros(1,NoIts);

theMSELMS = zeros(1,NoIts+1);

%=================================================================================
k = [-12.5 -2.5 -0.5 -0.2  0  0.2  0.5  2.5  12.5]; 
Lenc1 = length(c1);
c2Both = zeros(length(k),Lenc1);
for kInd = 1:1:length(k)
   c1Mat(kInd,1:Lenc1) = c1;
end

figure(864);
clf
subplot(321)
subplot(322)
subplot(323)
subplot(324)
subplot(325)
subplot(326)

NoPlots = 6;

% ====================================================================================

for MCtr = 1:1:NoIts; %=======START LOOP HERE!========================================
   
ErrLMS(MCtr) = (PC1*(DataVec(MCtr+1))+PC2*(DataVec(MCtr))-c1EstLMS(MCtr)*(DataVec(MCtr+1))-c2EstLMS(MCtr)*(DataVec(MCtr)));
theMSELMS(MCtr) = ErrLMS(MCtr)^2;
theMSELMSActual(MCtr) = (PC1 - c1EstLMS(MCtr))^2 + (PC2 - c2EstLMS(MCtr))^2;

%set(displayMSE,'string',['MSE = ',num2str(theMSELMSActual(MCtr),3),' at Iteration ',num2str(MCtr)],'fontsize',[10],'fontweight','demi')

PowerInFilt = DataVec(MCtr+1)^2 + DataVec(MCtr)^2 + 0.1;

c1EstLMS(1,MCtr+1) = c1EstLMS(1,MCtr) + 2*Mu*ErrLMS(MCtr)*DataVec(MCtr+1)/PowerInFilt;
c2EstLMS(1,MCtr+1) = c2EstLMS(1,MCtr) + 2*Mu*ErrLMS(MCtr)*DataVec(MCtr)/PowerInFilt;


 if ~(MCtr-(NoIts-NoPlots)<1)  
   if MCtr==NoIts-NoPlots+1
      subplot(321)
   elseif MCtr==NoIts-NoPlots+2
      subplot(322)
   elseif MCtr==NoIts-NoPlots+3
      subplot(323)
   elseif MCtr==NoIts-NoPlots+4
      subplot(324)
   elseif MCtr==NoIts-NoPlots+5
      subplot(325)
   elseif MCtr==NoIts-NoPlots+6
      subplot(326) 
   end
   
   hold on
% coefficient locus before updating,i.e, heading toward latest goal
plot(c1EstLMS(1,1:MCtr),c2EstLMS(1,1:MCtr),'o');
plot(c1EstLMS(1,1:MCtr),c2EstLMS(1,1:MCtr),':'); 
xlabel(['c1; Iteration = ',num2str(MCtr)])
ylabel(['c2'])
axis([PC1-c1Lim  PC1+c1Lim  PC2-c2Lim  PC2+c2Lim])
%====generate lines showing LQHP performance surface that forms next goal===============
TapRat = DataVec(MCtr+1)/DataVec(MCtr);

for kIndex = 1:1:length(k)
c2Both(kIndex,1:Lenc1) = (PC1 - c1)*TapRat + PC2 + sign(k(kIndex))*(abs(k(kIndex))^0.5)/DataVec(MCtr);
kMat(kIndex,1:Lenc1) = sign(k(kIndex))*k(kIndex)*ones(1,Lenc1);

if k(kIndex)==0
   line(c1, c2Both(kIndex,1:Lenc1));
 else
   line(c1, c2Both(kIndex,1:Lenc1));
end 
    plot(PC1,PC2,'b*');
    xlabel(['c1; Iteration = ',num2str(MCtr)])
    ylabel(['c2'])
	axis([PC1-c1Lim,  PC1+c1Lim,  PC2-c2Lim,  PC2+c2Lim])
end

%  replot  locus to include updated point
pause(1)

plot(c1EstLMS(1,1:MCtr+1),c2EstLMS(1,1:MCtr+1),'o');
plot(c1EstLMS(1,1:MCtr+1),c2EstLMS(1,1:MCtr+1)); 

end
hold off

end
%=======================================================================================








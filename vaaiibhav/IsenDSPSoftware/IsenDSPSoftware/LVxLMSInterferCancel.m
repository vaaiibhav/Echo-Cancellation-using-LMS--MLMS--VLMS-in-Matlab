function LVxLMSInterferCancel(k,Mu,DHorSH,MuteDesSig,DblTkUpD,NoSampsVOX,VOXThresh,GP2,GP4,DP1,DP2,DP3,DP4,NumTaps)
% function LVxLMSInterferCancel(k,Mu,DHorSH,MuteDesSig,DblTkUpD,NoSampsVOX,... 
% VOXThresh,GP2,GP4,DP1,DP2,DP3,DP4,NumTaps)
%
% k is the amount of noise to add to the interference signal ('whoknowsSR8K.wav')
% Mu is the usual LMS update term weight; 
% DHorSH yields a Dual-H system if passed as 1, or a Single-H system if passed as 0.  
% MuteDesSig mutes the desired signal ('drwatsonSR8K.wav') at its very beginning if passed as 1, 
% about one-third of the way through if passed as 2, or no muting at all if passed as 0;
% DblTkUpD, if passed as 0, allows coefficient updating any time, but if
% passed as 1 prevents coefficient update when the most recent NoSampsVOX samples of
% the desired signal level are above VOXThresh.
% GP2 is the gain in the path from the desired sound source to the noise
% reference microphone, and GP4 is the gain in the path from the noise source
% to the desired signal microphone; 
% DP1 is the Delay in samples from the desired sound source to the desired sound source microphone; 
% DP2 is the Delay from the desired sound source to the noise microphone; 
% DP3 is the Delay from the noise source to the noise source microphone;
% DP4 is the Delay from the noise source to the desired sound microphone.
% NumTaps is the number of Delays or taps to use in the adaptive filter.
% The final filtered output signal can be played using the call 
% sound(ActualErr,8000) after making the call global ActualErr in the Command window
% Typical calls might be:
% LVxLMSInterferCancel(0.02,0.3,1,0,0,50,0.03,1,1,1,6,1,6,5) % cancels interference only
% LVxLMSInterferCancel(0.02,0.3,1,0,0,50,0.03,0.16,0.16,7,42,7,42,35)
% LVxLMSInterferCancel(0.02,0.3,1,0,1,50,0.03,0.06,0.06,1,6,1,6,5)
% LVxLMSInterferCancel(0.02,0.3,1,0,0,50,0.03,1,1,6,1,1,6,5) % cancels desired sig also
% LVxLMSInterferCancel(0.02,0.3,1,1,0,50,0.03,1,1,1,6,1,6,5) % best: D-H w/immediate muting
% LVxLMSInterferCancel(0.02,0.2,0,2,1,80,0.03,1,1,1,6,1,6,5) % Single-H
% LVxLMSInterferCancel(0.02,0.3,0,1,1,50,0.03,1,1,1,6,1,6,5) % Single-H
%   
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

global ActualErr
ActualErr = [];

VOXThresh = abs(real(VOXThresh));
NoSampsVOX = abs(real(NoSampsVOX));
if NoSampsVOX > 1000
   NoSampsVOX = 1000;
end
DHorSH = abs(real(DHorSH));
if ~(DHorSH==1|DHorSH==0)
   DHorSH=1;
end
k = real(k);
Mu = abs(real(Mu));

MuteDesSig = real(MuteDesSig);

if ~(MuteDesSig==0|MuteDesSig==1|MuteDesSig==2)
   MuteDesSig = 1;
end

if ~(DblTkUpD==1|DblTkUpD==0)
DblTkUpD = 1;
end

GP2 = abs(real(GP2));
GP4 = abs(real(GP4));

DP1 = abs(real(DP1));
DP2 = abs(real(DP2));
DP3 = abs(real(DP3));
DP4 = abs(real(DP4));

DP1 = min([DP1 1000]);
DP2 = min([DP2 1000]);
DP3 = min([DP3 1000]);
DP4 = min([DP4 1000]);

NeededTaps = (DP4-DP3);

%================================================================================
[y,Fs,bits] = wavread('drwatsonSR8K.wav');
drwatSR = Fs

y = y';
lenFile = length(y);
DesiredSignal = y/(max(abs(y)));
lenDesiredSignal = length(DesiredSignal);
DesiredSignal = (DesiredSignal/(max(abs(DesiredSignal)))); 

if MuteDesSig==1 % mutes immediately
DesiredSignal(1,1:fix(0.05*lenDesiredSignal)) = 0.001*randn(1,fix(0.05*lenDesiredSignal));
elseif MuteDesSig==2 % waits before muting
   lenMute = fix(0.275*lenDesiredSignal) - fix(0.2*lenDesiredSignal);
   DesiredSignal(1,fix(0.32*lenDesiredSignal):fix(0.395*lenDesiredSignal)) = 0.001*randn(1,fix(lenMute)+1);
   end

[y1,Fs2,bits2] = wavread('whoknowsSR8K.wav'); 
NoiseVec = y1';
NoiseVec = NoiseVec(1,1:lenDesiredSignal);
a = length(NoiseVec);
NoiseVec = NoiseVec + k*randn(1,a);
NoiseVec = (NoiseVec/(max(abs(NoiseVec))));
% PlayDataVecIntCanc = NoiseVec;
OrigNoiseVec = NoiseVec;

NoiseVec(1,DP3+1:a) = OrigNoiseVec(1,1:a-DP3);
NoiseVec(1,1:DP3) = 0;
NoiseVec(1,DP2+1:a) = NoiseVec(1,DP2+1:a) + GP2*DesiredSignal(1,1:a-DP2);

DesiredSignal(1,DP1+1:a)= DesiredSignal(1,1:a-DP1);
DesiredSignal(1,1:DP1) = 0;
DesiredSignal(1,DP4+1:a) = DesiredSignal(1,DP4+1:a)+ GP4*OrigNoiseVec(1,1:a-DP4);
OrigNoiseVec = []; 
%================================================================================

BestTapWt = zeros(length(NoiseVec),NumTaps);
ActualFiltOut = zeros(1,a);
ActualErr = zeros(1,length(NoiseVec));
NewCoeff = zeros(1,a);

n = 1:1:NumTaps;
NoSampsVec = 1:1:NoSampsVOX;
Litn = 1:1:10;
FiltOut(1,length(NoiseVec)) = 0;
ActualErr(1,length(NoiseVec)) = 0;
TapToPlot = zeros(1,length(NoiseVec));
NewTestTapWt = zeros(1,NumTaps);
TestErr = zeros(1,length(NoiseVec));
TestErr(1,1:NumTaps+Litn(length(Litn))) = ones(1,NumTaps+Litn(length(Litn)));
NewCoeff(1,length(NoiseVec)) = 0;
CurBestMSE(1,length(NoiseVec)) = 0;
CurBestMSE(1,1:NumTaps+Litn(length(Litn))) = ones(1,NumTaps+Litn(length(Litn)));
BestTapWt(1:length(NoiseVec),1:NumTaps) = 0;
TestERLE = 0.001;

for CurDtaPtr = NumTaps+Litn(length(Litn))+1:1:length(NoiseVec); 
   if DHorSH==0  % Single-H 
     ActualFiltOut =  sum(BestTapWt(CurDtaPtr,1:NumTaps).*NoiseVec(1,CurDtaPtr-n));
	 ActualErr(1,CurDtaPtr) =  DesiredSignal(1,CurDtaPtr) - ActualFiltOut;   
   %========================================================================
 if DblTkUpD==1
   if CurDtaPtr > NoSampsVOX    
      if sqrt((1/NoSampsVOX)*sum(DesiredSignal(1,CurDtaPtr-NoSampsVec).^2))< VOXThresh 
         BestTapWt(CurDtaPtr+1,1:NumTaps) = BestTapWt(CurDtaPtr,1:NumTaps) + ...
   		2*(Mu)*ActualErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n))/sum((NoiseVec(1,CurDtaPtr-n)).^2); 
      NewCoeff(1,CurDtaPtr) = 1;
   	else
      BestTapWt(CurDtaPtr+1,1:NumTaps) = BestTapWt(CurDtaPtr,1:NumTaps);
      end
   else
      BestTapWt(CurDtaPtr+1,1:NumTaps) = BestTapWt(CurDtaPtr,1:NumTaps) + ...
   		2*(Mu)*ActualErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n))/sum((NoiseVec(1,CurDtaPtr-n)).^2); 
      NewCoeff(1,CurDtaPtr) = 1;
   end
else % Update all the time even if doubletalk
BestTapWt(CurDtaPtr+1,1:NumTaps) = BestTapWt(CurDtaPtr,1:NumTaps) + ...
   2*(Mu)*ActualErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n))/sum((NoiseVec(1,CurDtaPtr-n)).^2);
NewCoeff(1,CurDtaPtr) = 1;
end

CurBestMSE(1,CurDtaPtr) = sum(NoiseVec(1,CurDtaPtr-Litn).^2)/(sum(ActualErr(1,CurDtaPtr-Litn).^2)+10^-6);   
    
else   %  DUAL-h STARTS HERE    
         if TestERLE > CurBestMSE(1,CurDtaPtr-1)
            BestTapWt(CurDtaPtr,1:NumTaps) = NewTestTapWt;
      	CurBestMSE(1,CurDtaPtr) = TestERLE;
      	else
            CurBestMSE(1,CurDtaPtr) = CurBestMSE(1,CurDtaPtr-1);
            BestTapWt(CurDtaPtr,1:NumTaps) = BestTapWt(CurDtaPtr-1,1:NumTaps);
   		end  
ActualFiltOut =  sum(BestTapWt(CurDtaPtr,1:NumTaps).*NoiseVec(1,CurDtaPtr-n));
TestFiltOut =  sum(NewTestTapWt.*NoiseVec(1,CurDtaPtr-n));
ActualErr(1,CurDtaPtr) =  DesiredSignal(1,CurDtaPtr) - ActualFiltOut;          
TestErr(1,CurDtaPtr) =  DesiredSignal(1,CurDtaPtr) - TestFiltOut;   
TestERLE = sum(NoiseVec(1,CurDtaPtr-Litn).^2)/(sum(TestErr(1,CurDtaPtr-Litn).^2)+10^-6);   
if DblTkUpD==1
   if CurDtaPtr> NoSampsVOX   
       if sqrt((1/NoSampsVOX)*sum(DesiredSignal(1,CurDtaPtr-NoSampsVec).^2))< VOXThresh
         NewTestTapWt = NewTestTapWt + ...
   		2*Mu*TestErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n))/sum((NoiseVec(1,CurDtaPtr-n)).^2); 
			NewCoeff(1,CurDtaPtr) = 1;
      end
   else
    		NewTestTapWt = NewTestTapWt + ...
   		2*Mu*TestErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n))/sum((NoiseVec(1,CurDtaPtr-n)).^2); 
			NewCoeff(1,CurDtaPtr) = 1;  
   end
else % Update all the time even if doubletalk
	NewTestTapWt = NewTestTapWt + ...
   2*Mu*TestErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n))/sum((NoiseVec(1,CurDtaPtr-n)).^2);
	NewCoeff(1,CurDtaPtr) = 1;
end

end
end
%==========================================================================================

figure(45)
clf
subplot(411)

line([0 a],[1 1]);
  hold on
line([0 a],[0 0]);

   for ctr = 1:1:NumTaps
      plot(BestTapWt(:,ctr))
	end
  hold off
xlabel(['(a) Sample Number'])
ylabel(['Amplitude'])
plotlim = 1.2*max(max(BestTapWt))+ 0.05;
plotlimlow = min(min(BestTapWt))-0.05;
axis([0 length(NoiseVec) plotlimlow  plotlim])

subplot(412)
plot(ActualErr(1,1:length(NoiseVec)))
plotlim1 = max(abs(ActualErr));
axis([0 length(NoiseVec) -1.2*plotlim1 1.2*plotlim1])
ylabel(['Amplitude'])
xlabel(['(b) Sample, Error/Output'])

subplot(413)
hndPlot = plot(NewCoeff(1,:),'b'); 
ylabel(['Binary'])
xlabel(['(c) Sample, Fcn Coefficient Update Permitted?'])
axis([0 length(NewCoeff) -0.2  1.2]) 

subplot(414)
plot(20*log10(CurBestMSE(1,1:length(CurBestMSE))+10^-6))
plotmin = 1.2*min(20*log10( CurBestMSE(1,1:length(CurBestMSE))+10^-6));
plotmax = 1.2*max(20*log10(CurBestMSE(1,1:length(CurBestMSE))+10^-6));
xlabel(['(d) Sample, CurBestIRLE'])
ylabel(['Mag, dB'])
axis([0 length(CurBestMSE) plotmin  plotmax])

%==================================================================
ActualErr = (0.99/max(abs(ActualErr)))*ActualErr;
sound(ActualErr,Fs)
pause(5)
sound(ActualErr,Fs) % repeat





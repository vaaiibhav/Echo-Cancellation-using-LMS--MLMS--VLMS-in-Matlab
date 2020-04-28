function LVxLMSInterferCancShort(k,Mu,DHorSH,MuteDesSig,GP2,GP4,DP1,DP2,DP3,DP4,NumTaps,PrPC)
% function LVxLMSInterferCancShort(k,Mu,DHorSH,MuteDesSig,GP2,GP4,DP1,DP2,DP3,DP4,NumTaps,PrPC)
% k is the amount of noise to add to the interference signal ('whoknowsSR8K.wav')
% Mu is the usual LMS update term weight; 
% DHorSH yields a Dual-H system if passed as 1, or a Single-H system if passed as 0.  
% MuteDesSig mutes the desired signal ('drwatsonSR8K.wav') at its very beginning if passed as 1, 
% about one-third of the way through if passed as 2, or no muting at all if
% passed as 0;
% GP2 is the gain in the path from the desired sound source to the noise
% reference microphone, and GP4 is the gain in the path from the noise source
% to the desired signal microphone; 
% DP1 is the Delay in samples from the desired sound source to the desired sound source microphone; 
% DP2 is the Delay from the desired sound source to the noise microphone; 
% DP3 is the Delay from the noise source to the noise source microphone;
% DP4 is the Delay from the noise source to the desired sound microphone.
% NumTaps is the number of Delays or taps to use in the adaptive filter.
% The total amount of computation can be limited by processing only a
% percentage of the audio signals, and that percentage is passed as the
% input variable PrPC.
% whoknowsSR4K.wav and drwatsonSR4K.wav as test signals, 100%=22184 samps. 
% The final filtered output signal is interpolated by a factor of 2 using the 
% function interp to raise its sample rate to 8000 Hz so it can be played using the
% call sound(ActualErr,8000) after making the call global ActualErr in the
% Command window
%
% Typical calls might be:
% LVxLMSInterferCancShort(0.02,0.3,1,0,1,1,1,6,1,6,5,30) % cancels interference only
% LVxLMSInterferCancShort(0.02,0.3,1,0,0.16,0.16,7,42,7,42,35,30)
% LVxLMSInterferCancShort(0.02,0.3,1,0,0.06,0.06,1,6,1,6,5,30)
% LVxLMSInterferCancShort(0.02,0.3,1,0,1,1,6,1,1,6,5,30)  % cancels desired sig also
% LVxLMSInterferCancShort(0.02,0.3,1,1,1,1,1,6,1,6,5,30)  % best: D-H w/immediate muting
% LVxLMSInterferCancShort(0.02,0.2,0,2,1,1,1,6,1,6,5,30)  % Single-H
% LVxLMSInterferCancShort(0.02,0.3,0,1,1,1,1,6,1,6,5,30)  % Single-H
%   
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

global ActualErr
global DesiredSignal

ActualErr = [];
DesiredSignal = [];

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

ActualErr = [];

%================================================================================
[y,Fs,bits] = wavread('drwatsonSR4K.wav');

y = y';
lenFile = length(y);
len2Do = fix(lenFile*PrPC/100);
DesiredSignal = y(1:len2Do)/(max(abs(y(1:len2Do))));
lenDesiredSignal = length(DesiredSignal);
DesiredSignal = (DesiredSignal/(max(abs(DesiredSignal)))); 

if MuteDesSig==1 % mutes immediately
DesiredSignal(1,1:fix(0.05*lenDesiredSignal)) = 0.001*randn(1,fix(0.05*lenDesiredSignal));
elseif MuteDesSig==2 % waits before muting
   lenMute = fix(0.275*lenDesiredSignal) - fix(0.2*lenDesiredSignal);
   DesiredSignal(1,fix(0.32*lenDesiredSignal):fix(0.395*lenDesiredSignal)) = 0.001*randn(1,fix(lenMute)+1);
   end

[y1,Fs2,bits2] = wavread('whoknowsSR4K.wav'); 
NoiseVec = y1';
NoiseVec = NoiseVec(1,1:lenDesiredSignal);
lenNV = length(NoiseVec);
NoiseVec = NoiseVec + k*randn(1,lenNV);
NoiseVec = (NoiseVec/(max(abs(NoiseVec))));
OrigNoiseVec = NoiseVec;

NoiseVec(1,DP3+1:lenNV) = OrigNoiseVec(1,1:lenNV-DP3);
NoiseVec(1,1:DP3) = 0;
NoiseVec(1,DP2+1:lenNV) = NoiseVec(1,DP2+1:lenNV) + GP2*DesiredSignal(1,1:lenNV-DP2);

DesiredSignal(1,DP1+1:lenNV)= DesiredSignal(1,1:lenNV-DP1);
DesiredSignal(1,1:DP1) = 0;
DesiredSignal(1,DP4+1:lenNV) = DesiredSignal(1,DP4+1:lenNV)+ GP4*OrigNoiseVec(1,1:lenNV-DP4);
OrigNoiseVec = []; 
%==========================================================================
n = 1:1:NumTaps;
Litn = 1:1:10;
ActualFiltOut = zeros(1,lenNV);
ActualErr = zeros(1,length(NoiseVec));
NewTestTapWt = zeros(1,NumTaps);
TestErr = zeros(1,length(NoiseVec));
TestErr(1,1:NumTaps+Litn(length(Litn))) = ones(1,NumTaps+Litn(length(Litn)));
CurBestMSE = 1;
BestTapWt = zeros(1,NumTaps);
TestERLE = 0.001;

% compute update term normalizing weight efficiently
sqNV = NoiseVec.^2;
normMag = filter(ones(1,max(n)),1,sqNV);
upDateWt = 2*Mu./normMag;

for CurDtaPtr = NumTaps+Litn(length(Litn))+1:1:length(NoiseVec); 
    theCtr = CurDtaPtr
  if DHorSH==0  % Single-H 
    ActualFiltOut =  sum(BestTapWt.*NoiseVec(1,CurDtaPtr-n));
	ActualErr(1,CurDtaPtr) =  DesiredSignal(1,CurDtaPtr) - ActualFiltOut;   
    BestTapWt = BestTapWt + ...
    upDateWt(CurDtaPtr)*ActualErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n));
CurBestMSE = normMag(CurDtaPtr)/(sum(ActualErr(1,CurDtaPtr-Litn).^2)+10^-6);      
else   %  DUAL-h STARTS HERE    
    if TestERLE > CurBestMSE
       BestTapWt = NewTestTapWt;
      CurBestMSE = TestERLE;
%   else
   	end  
ActualFiltOut =  sum(BestTapWt.*NoiseVec(1,CurDtaPtr-n));
TestFiltOut =  sum(NewTestTapWt.*NoiseVec(1,CurDtaPtr-n));
ActualErr(1,CurDtaPtr) =  DesiredSignal(1,CurDtaPtr) - ActualFiltOut;          
TestErr(1,CurDtaPtr) =  DesiredSignal(1,CurDtaPtr) - TestFiltOut;   
TestERLE = normMag(CurDtaPtr)/(sum(TestErr(1,CurDtaPtr-Litn).^2)+10^-6);   
NewTestTapWt = NewTestTapWt + ...
upDateWt(CurDtaPtr)*TestErr(1,CurDtaPtr)*(NoiseVec(1,CurDtaPtr-n));
end
end
%==========================================================================================

figure(45)
subplot(211)
plot(DesiredSignal,'k')
plotlim1 = max(abs(DesiredSignal));
ylabel(['Amplitude'])
xlabel(['(a) Desired Signal, Plus Noise'])
axis([0, length(DesiredSignal), -1.2*plotlim1, 1.2*plotlim1])

subplot(212)
plot(ActualErr,'k')
plotlim1 = max(abs(ActualErr));
axis([0, length(ActualErr), -1.2*plotlim1, 1.2*plotlim1])
ylabel(['Amplitude'])
xlabel(['(b) Error/Output Sample Number'])

%==================================================================
DesiredSignal = interp(DesiredSignal,2,10,0.5);
DesiredSignal = (0.99/max(abs(DesiredSignal)))*DesiredSignal;
ActualErr = interp(ActualErr,2,10,0.5);
ActualErr = (0.99/max(abs(ActualErr)))*ActualErr;
sound(DesiredSignal,2*Fs)
pause(5)
sound(ActualErr,2*Fs) 

Comment = 'global output sound variables are DesiredSignal and ActualErr'
Comment = ['Output sound variable sample rate = ',num2str(2*Fs)]





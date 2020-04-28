function LVxLMSAdptFiltEchoShort(k,Mu,DHorSH,MuteNrEnd,PrPC)
% function LVxLMSAdptFiltEchoShort(k,Mu,DHorSH,MuteNrEnd,TstSig,PrPC)
% k is the amount of noise to add to the Far End signal (see below)
% Mu is the usual LMS update term weight;
% DHorSH yields a Dual-H system if passed as 1, or a Single-H system if passed as 0; 
% MuteNrEnd mutes the desired signal ('drwatsonSR4K') at its very beginning
% if passed as 1, about one-third of the way through if passed as 2, or no muting at all if
% passed as 0. 
% The Near End signal consists of the audio file 'drwatsonSR4K.wav'
% The Far End signal consists of the audio file 'whoknowsSR4K.wav' plus 
% random noise weighted by k, limited to the length of the Near End signal;
% The total amount of computation can be limited by processing only a
% percentage of the audio signals, and that percentage is passed as the
% input variable PrPC.
% whoknowsSR4K.wav and drwatsonSR4K.wav as test signals, 100%=22184 samps
% The adaptive FIR is ten samples long and the echo is simulated as a single delay of 6 samples, for 
% example. The final filtered output signal is interpolated by a factor of 2 using the function interp
% to raise its sample rate to 8000 Hz so it may be played using the call sound(EchoErr,8000) after 
% making the call global EchoErr in the Command window.
% Test calls:
% LVxLMSAdptFiltEchoShort(0.1,0,0,0,50) % Mu = 0, Far End easily heard
% LVxLMSAdptFiltEchoShort(0.1,0.2,0,0,50) % S-H, no Mute
% LVxLMSAdptFiltEchoShort(0.1,0.2,0,1,50) % S-H, Mute Imm
% LVxLMSAdptFiltEchoShort(0.1,0.2,1,0,50) % D-H, no Mute
% LVxLMSAdptFiltEchoShort(0.1,0.2,1,1,50) % D-H, Mute Imm
% LVxLMSAdptFiltEchoShort(0.1,0.2,1,2,50) % D-H, Mute delayed
% LVxLMSAdptFiltEchoShort(0.01,0.2,0,0,50) % S-H, no Mute
% LVxLMSAdptFiltEchoShort(0.01,0.2,0,1,50) % S-H, Mute Imm
% LVxLMSAdptFiltEchoShort(0.01,0.2,1,0,50) % D-H, no Mute
% LVxLMSAdptFiltEchoShort(0.01,0.2,1,1,50) % D-H, Mute Imm
% LVxLMSAdptFiltEchoShort(0.01,0.2,1,2,50) % D-H, Mute delayed
% LVxLMSAdptFiltEchoShort(0.01,0.02,1,0,100) % Low Mu, can hear Far End
% gradually diminish
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
global EchoErr
EchoErr = [];

DHorSH = abs(real(DHorSH));
if ~(DHorSH==1|DHorSH==0)
   DHorSH=1;
end
k = real(k);
Mu = abs(Mu);

MuteNrEnd = real(MuteNrEnd);

if ~(MuteNrEnd==0|MuteNrEnd==1|MuteNrEnd==2)
   MuteNrEnd = 1;
end

InterpFac = 2;
[NrEndSpeech,Fs,bits] = wavread('drwatsonSR4K.wav');
lendrwatson = length(NrEndSpeech);
Len2Do = fix(lendrwatson*PrPC/100);
NrEndSpeech = NrEndSpeech(1:Len2Do);
NrEndSpeech = NrEndSpeech/(max(abs(NrEndSpeech)));
lenNearEnd = length(NrEndSpeech);

if MuteNrEnd==1 % mutes immediately
NrEndSpeech(1:fix(0.05*lenNearEnd),1) = 0.001*randn(fix(0.05*lenNearEnd),1);
elseif MuteNrEnd==2 % waits before muting
 lenMute = fix(0.25*lenNearEnd) - fix(0.2*lenNearEnd);
 NrEndSpeech(fix(0.32*lenNearEnd):fix(0.37*lenNearEnd),1) = 0.001*randn(lenMute+1,1);
end

[FarEnd,Fs2,bits2] = wavread('whoknowsSR4K.wav');
FarEnd = FarEnd(1:Len2Do);
FarEnd = FarEnd(1:lenNearEnd,1);
lenFarEnd = length(FarEnd);
FarEnd  = FarEnd  + k*randn(lenFarEnd,1);
FarEnd = FarEnd/max(abs(FarEnd));

%==========================================================
Litn = 1:1:10;
TapNo = 6;
FiltOut = 0;
EchoErr = zeros(lenFarEnd,1);
n = 1:1:10; 
BestTapWt = zeros(10,1);
NewTestTapWt = zeros(10,1);
TestErr = zeros(lenFarEnd,1);
TestErr = ones(11+TapNo+Litn(length(Litn)),1);
CurBestMSE = 1; 
TestERLE = 0.00001;

dtaPtr = TapNo:1:lenNearEnd;
NetNrEndPlFar = FarEnd(dtaPtr-TapNo+1,1) + NrEndSpeech(dtaPtr,1);

% compute update term normalizing weight efficiently
sqFE = FarEnd.^2;
normMag = filter(ones(1,max(n)),1,sqFE);
upDateWt = 2*Mu./normMag;
%===================================================
LoLim = max([11+TapNo+Litn(length(Litn))+1]);
 
if DHorSH==0  % Single-H 
    for CurDtaPtr = LoLim:1:lenNearEnd-10;   
ccc = FarEnd(CurDtaPtr+1-n,1);
FiltOut = sum(BestTapWt.*ccc);
ddd = NetNrEndPlFar(CurDtaPtr,1) - FiltOut; 
EchoErr(CurDtaPtr,1) = ddd;
BestTapWt = BestTapWt + upDateWt(CurDtaPtr)*ddd*ccc;
CurBestMSE = normMag(CurDtaPtr)/(sum(ddd.^2)+10^(-6));  
    end
end

if DHorSH==1  % Dual-H  % =============Dual-H code here============    
CurDtaPtr = LoLim;
 while CurDtaPtr < lenNearEnd-10; 
        if TestERLE > CurBestMSE
            BestTapWt = NewTestTapWt;
            CurBestMSE = TestERLE;
        end  
aaa = FarEnd(CurDtaPtr-n+1,1);
FiltOut =  sum(BestTapWt.*aaa);
EchoErr(CurDtaPtr,1) = NetNrEndPlFar(CurDtaPtr,1)-FiltOut;
TestFiltOut =  sum(NewTestTapWt.*aaa);
TestErr(CurDtaPtr,1) = NetNrEndPlFar(CurDtaPtr,1)-TestFiltOut;  
TestERLE = normMag(CurDtaPtr)/(sum(TestErr(CurDtaPtr-Litn+1,1).^2)+10^-6);    
NewTestTapWt = NewTestTapWt + upDateWt(CurDtaPtr)*TestErr(CurDtaPtr,1)*aaa;
CurDtaPtr = CurDtaPtr + 1;
end
end

figure(37)
subplot(211)
plot(NetNrEndPlFar)
plotlim2 = max(abs(NetNrEndPlFar));
xlabel(['(a) Sample, Near End Plus Echo'])
ylabel(['Amplitude'])
axis([0, lenFarEnd, -1.2*plotlim2, 1.2*plotlim2])

subplot(212)
plot(EchoErr)
plotlim1 = max(abs(EchoErr));
ylabel(['Amplitude'])
xlabel(['(b) Sample, Error/Output'])
axis([0, lenFarEnd, -1.2*plotlim1, 1.2*plotlim1])

EchoErr = (0.99/max(abs(EchoErr)))*EchoErr;

FinalSampleRate = Fs*InterpFac; 
EchoErr = interp(EchoErr,InterpFac,10,(1/InterpFac));
sound(EchoErr,FinalSampleRate)
Comment = 'global sound output variable is EchoErr & sample rate is 8000'




  






function LVxLMS2TapFIR(PC1,PC2, NoIts,tstDType,CosFrq,TwoMu,NoiseAmp)
% PC1 and PC2 are Plant Coefficients, to be modeled by the 2-tap FIR
% NoIts is the test sequence length in samples
% tstDType selects test signals to pass through the Plant and the filter,
% as follows:
% 1 = Random Noise; 2 = Unit Step; 3 = Nyquist Limit (Fs/2); 4 = Half-Band (Fs/4)
% 5 = cosine of freq CosFrq; 6 = Cos(1Hz) + Cos(Fs/4) + Cos(Fs/2)
% TwoMu is the update term weight
% Random noise having amplitude NoiseAmp is added to the test signal
% selected by tstDType.
% Test calls:
% LVxLMS2TapFIR(4,-3,64,1,[],0.5,0)
% LVxLMS2TapFIR(4,-3,64,2,[],0.5,0)
% LVxLMS2TapFIR(4,-3,64,3,[],0.5,0)
% LVxLMS2TapFIR(4,-3,64,4,[],0.5,0)
% LVxLMS2TapFIR(4,-3,64,5,[25],0.5,0)
% LVxLMS2TapFIR(4,-3,64,5,[1],0.5,0)
% LVxLMS2TapFIR(4,-3,64,6,[],0.5,0)
% LVxLMS2TapFIR(4,-3,64,0,[],1.15,0) 
% LVxLMS2TapFIR(4,-3,7,4,[],1.15,0)  
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
figure(17);
clf
seqLen = NoIts + 1;

if tstDType==1
 DataVec = randn(1,seqLen);
 SigType = 'Random Noise'; 
elseif tstDType==2
 DataVec = ones(1,seqLen);  %  DC
 SigType = 'Unit Step (DC)';
elseif tstDType==3
 DataVec = zeros(1,seqLen);
 nn = 1:1:seqLen/2;
 DataVec(1,2*nn ) = -1;
 DataVec(1,2*nn -1) = 1; % Nyquist limit frequency (fs/2)
 SigType = 'Nyquist Limit (fs/2)';
elseif tstDType==4
 DataVec = zeros(1,seqLen);
 nn = 1:1:seqLen/4;
 DataVec(1,4*nn ) = -1;
 DataVec(1,4*nn -2) = 1; % half-band frequency (fs/4)
 DataVec = DataVec +  0.001*randn(1,seqLen);
 SigType = 'Half-Band Freq (fs/4)';
elseif tstDType==5   % cosine      
DataVec = cos(2*pi*CosFrq*(1:1:seqLen)/seqLen); 
SigType = 'Cosine Wave';
elseif tstDType==6   % 3 cosines      
DataVec = cos(2*pi*1*(1:1:seqLen)/seqLen) + cos(2*pi*(seqLen/4)*(1:1:seqLen)/seqLen) + cos(2*pi*(seqLen/2)*(1:1:seqLen)/seqLen);
SigType = 'Three Cosine Waves';
else
 DataVec = randn(1,seqLen);
 SigType = 'Random Noise';
end

DataVec = DataVec/max(abs(DataVec));
DataVec = DataVec + NoiseAmp*randn(1,seqLen);

plotlim = 1.2*max(abs(DataVec));

subplot(221)
hold on
stem(DataVec,'bo');
ylabel(['Amplitude'])
xlabel(['(a) Input Sample'])
grid on
axis([0  NoIts  -plotlim  plotlim])

%=================================================================================
PlantCoeff1 = PC1;
PlantCoeff2 = PC2;
c1EstLMS = zeros(1,NoIts);
c2EstLMS = zeros(1,NoIts);
c1EstLMS(1,1) = 0;
c2EstLMS(1,1) = 0;
ErrLMS = zeros(1,NoIts);

for MCtr = 1:1:NoIts
ErrLMS(1,MCtr) = (PlantCoeff1*(DataVec(MCtr+1))+PlantCoeff2*(DataVec(MCtr))-c1EstLMS(MCtr)*(DataVec(MCtr+1))-c2EstLMS(MCtr)*(DataVec(MCtr)));
PowerInFilt = DataVec(MCtr+1)^2 + DataVec(MCtr)^2 + 0.1;
c1EstLMS(1,MCtr+1) = c1EstLMS(1,MCtr) + TwoMu*ErrLMS(MCtr)*DataVec(MCtr+1)/PowerInFilt;
c2EstLMS(1,MCtr+1) = c2EstLMS(1,MCtr) + TwoMu*ErrLMS(MCtr)*DataVec(MCtr)/PowerInFilt;
%=======================================================================================
if MCtr==NoIts
   break
end
end

ErrPlotLim = 1.2*max(abs(ErrLMS));
final_c1Est = c1EstLMS(1,length(c1EstLMS))
final_c2Est = c2EstLMS(1,length(c2EstLMS))

subplot(222)
if sign(ErrLMS).*abs(ErrLMS)==ErrLMS  % ErrLMS is real only
    theSig = ErrLMS(1,1:MCtr);
    xvec = 1:1:length(theSig);
stem(xvec,theSig,'bo');
xlabel(['(b) Sample/Iteration'])
ylabel(['Error'])
grid on
else
theSig = abs(ErrLMS(1,1:MCtr));
xvec = 1:1:length(theSig);
stem(xvec,theSig,'bo');   % ErrLMS complex
xlabel(['(b) Sample/Iteration'])
ylabel(['Error'])
grid on
end
axis([0 NoIts -ErrPlotLim ErrPlotLim])

subplot(223)

stem(c1EstLMS)
xlabel(['(c) Iteration'])
ylabel(['Amplitude, c1'])
grid on
axis([0, NoIts, -1.2*max(abs(c1EstLMS)), 1.2*max(abs(c1EstLMS))])

subplot(224)

stem(c2EstLMS)
xlabel(['(d) Iteration'])
ylabel(['Amplitude, c2'])
grid on
axis([0, NoIts, -1.2*max(abs(c2EstLMS)), 1.2*max(abs(c2EstLMS))])





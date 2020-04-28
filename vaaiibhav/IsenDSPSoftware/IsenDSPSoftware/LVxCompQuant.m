function LVxCompQuant(TypeCompr,AMu,NoBits)
%
% LVxCompQuant(1,255,3) % no comp
% LVxCompQuant(2,255,3) % mu-law compr
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

strDrReadPath = ['drwatsonSR8K.wav']; 
[AudioSig1,Fs,bits] = wavread(strDrReadPath); 
MuLawSR = 8000;
AudioSig1 = AudioSig1';
%AudioSig1 = resample(AudioSig1,8,22);

%AudioSig1 = AudioSig1 - mean(AudioSig1);
LenAudioSig1 = length(AudioSig1)
AudioSig1 = AudioSig1/(max(abs(AudioSig1)));
dur = (length(AudioSig1)/MuLawSR);

astart = 1;
aend = 795;
t = (0:1:fix(LenAudioSig1/5))/fix(LenAudioSig1)*dur;
LenOft = length(t)

theRawInputSound = AudioSig1;
strTestSig = 'Speech & Music';

scalfac = (max(abs(theRawInputSound)));
theRawInputSound = theRawInputSound/scalfac;
theRawInputSound = 1*theRawInputSound; 

figure(17)
clf

subplot(311)
plot(theRawInputSound);
xlabel(['Test Signal: ',strTestSig])
ylabel(['Amplitude'])
axis([0  length(theRawInputSound)  min(theRawInputSound)-0.1   max(theRawInputSound)+0.1])

Quantlev = NoBits;
 
zx = 1:1:length(theRawInputSound);

TypeAMu = TypeCompr;
if TypeAMu==1
strTypeComp = 'NoCompr';
  % no compression
theRawInputSoundcomp = theRawInputSound;
theRawInputSoundcomp = theRawInputSoundcomp/(max(abs(theRawInputSoundcomp)));
Gran = 0.9999*(1/((2^Quantlev)-1)); 
qtheRawInputSoundcomp(zx) = (fix(theRawInputSoundcomp(zx)/Gran))*Gran;
Decomp = qtheRawInputSoundcomp;

subplot(312)
plot(Decomp)   
xlabel(['Uncompressed, Quantized to ',num2str(Quantlev),' bits'])
 plotlim = 1.1*max(abs(Decomp));
 ylabel(['Amplitude'])
 axis([0  length(Decomp)   (min(Decomp)-0.1)  max([plotlim  1.1])    ])

subplot(313)
plot(Decomp)   
xlabel(['Output/Decompressed Signal']) 
plotlim = 1.1*max(abs(Decomp));
ylabel(['Amplitude'])
 axis([0 length(Decomp)  (min(Decomp)-0.1)  max([plotlim  1.1]) ])

theCompQuantPlaySig = Decomp;

elseif TypeAMu==2; % Mu-law compression
strTypeComp = 'MuLaw';
          
theRawInputSoundcomp(zx)=((log(1+AMu*(abs(theRawInputSound(zx)))))/(log(1+AMu))).*sign(theRawInputSound(zx));
theRawInputSoundcomp = theRawInputSoundcomp/(max(abs(theRawInputSoundcomp)));        
Gran = 0.99*(1/((2^Quantlev)-1)); 
 
qtheRawInputSoundcomp(zx) = (fix(theRawInputSoundcomp(zx)/Gran))*Gran;

C1 = log(1+AMu);
% decompress (expand) the quantized version
theexp = (abs(qtheRawInputSoundcomp(zx)));
partans = exp(C1*theexp);
Decomp(zx)=sign(qtheRawInputSoundcomp(zx)).*(partans -1)/AMu;
% decompress the compressed analog signal
Antheexp = (abs(theRawInputSoundcomp(zx)));
Anpartans = exp(C1*Antheexp);
AnalogDecomp(zx)=sign(theRawInputSoundcomp(zx)).*(Anpartans -1)/AMu;
%==================================================================
subplot(312)
plot(qtheRawInputSoundcomp)   
xlabel(['Compressed by Mu-Law, Quantized to ',num2str(Quantlev),' bits'])
plotlim = 1.1*max(abs(qtheRawInputSoundcomp));
ylabel(['Amplitude'])
axis([0  length(qtheRawInputSoundcomp)   (min(qtheRawInputSoundcomp)-0.1)  max([plotlim  1.1])    ])
 
subplot(313)
plot(Decomp)   
xlabel(['Decompressed (Expanded) Signal']) 
plotlim = 1.1*max(abs(Decomp));
ylabel(['Amplitude'])
axis([0 length(Decomp)  (min(Decomp)-0.1)  max([plotlim  1.1]) ])
else
    error('Pass TypeCompr as 1 for no compression or 2 for Mu-law compression')
end

sound(Decomp,MuLawSR)




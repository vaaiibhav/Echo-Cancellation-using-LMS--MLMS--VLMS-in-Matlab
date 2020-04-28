function LVPaddedDFTMovie(ImpResp,LenDTFT, theCompMode)
% LVPaddedDFTMovie([1 1 1 1],128,1)
% LVPaddedDFTMovie([1 0 0 1],128,1)
% LVPaddedDFTMovie([1 0 0 0 0 0 0 1],128,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

if max(abs(ImpResp))>1.25
   ImpResp = (1.25/max(abs(ImpResp)))*ImpResp;
end

SR = LenDTFT;
ImpResp = [ImpResp zeros(1,(SR-length(ImpResp)))];

Mag = zeros(1,fix(SR/2));

t = 0:1/SR:1-1/SR;

figure(640)
clf

for TestFreq = 0:1:fix(SR/2);
   
TestFreqReal = cos(2*pi*TestFreq*t);

TestFreqImag = -sin(2*pi*TestFreq*t);

subplot(311)

impToPlot = find(ImpResp);
stem(impToPlot-1,ImpResp(impToPlot),'b*');
hold on;
stem(0:1:length(TestFreqReal)-1,TestFreqReal,'ro');
plot(0:1:length(TestFreqReal)-1,TestFreqReal,'r')
axis([0 LenDTFT -1.3 2.75])
hold off;

xlabel(['(a)  Sample'])
ylabel(['Amp'])

ComplexCorr = (1/SR)*sum((TestFreqReal +j*TestFreqImag).*ImpResp);
realCorrVal = real(ComplexCorr);

if abs(realCorrVal)<10^-10
   realCorrVal = 0;
end
%text(5,2,['Correlation Value = ',num2str(realCorrVal)])

subplot(312)

stem(impToPlot-1,ImpResp(impToPlot),'b*');
hold on;
stem(0:1:length(TestFreqImag)-1,TestFreqImag,'ro');
plot(0:1:length(TestFreqImag)-1,TestFreqImag,'r')
axis([0 LenDTFT -1.3 2.75])
hold off;
xlabel(['(b)  Sample'])
ylabel(['Amp'])

imagCorrVal = imag(ComplexCorr);
if abs(imagCorrVal)<10^-10
   imagCorrVal = 0;
end

%text(5,2,['Correlation Value = ',num2str(imagCorrVal)])

if TestFreq==0|TestFreq==fix(SR/2)
Mag(1,TestFreq+1) = sqrt( (realCorrVal)^2 + (imagCorrVal)^2 );
else
Mag(1,TestFreq+1) = sqrt( (realCorrVal)^2 + (imagCorrVal)^2 );
end

MaxMag = 1.2*max(abs(fft(ImpResp,LenDTFT)))/LenDTFT;

subplot(313)
stem(0:1:length(Mag)-1,Mag,'bo');
xlabel(['(c)  Bin'])
ylabel(['Mag'])
axis([0 fix(LenDTFT/2) 0 MaxMag])

if TestFreq == fix(SR/2)
   return
end

if theCompMode==1
pause(0.2)
else
    pause
end
end
function BinOut = LVxBinaryCodeMethods(BitsQ,SR,Bias,Freq,Amp,CodeMeth,PlotType)
% Quantizes a sine wave of amplitude Amp and frequency Freq using BitsQ
% number of quantization bits at a sample rate of SR. 
% Bias: 0 for none, 1 for 1/2 LSB
% CodeMeth: 1 % Sign + Mag; CodeMeth==2 = Offset;
% PlotType: val 0 plot as multiples of LSB, or 1 to plot in volts
% Generates a figure with three subplots, the first is the (simulated) analog test signal, the second
% has an overlay of the simulated test signal and its quantized version, and the third is the quantization
% error
% Test calls:
% BinOut = LVxBinaryCodeMethods(4, 1000,0, 10,170,2,0)
% BinOut = LVxBinaryCodeMethods(2, 1000,0, 10,170,1,1)
% BinOut = LVxBinaryCodeMethods(3, 1000,1, 10,170,1,1)
% BinOut = LVxBinaryCodeMethods(3, 1000,1, 10,170,2,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
UpperPlotAnalogYes1 = 1;

figure(9586);
clf

%========simplesampling compute=============
valNormalOrNo = PlotType;
UpperPlotAnalogYes1 = 1;
        
t = 0:1/SR:(1/10 - 1/SR);

y = Amp*sin(2*pi*t*Freq);

subplot(311)
if UpperPlotAnalogYes1==0
   stem(t,y,'bo')
   xlabel(['(a)  A sampled ',num2str(Freq),' Hz sine wave.  Sample rate ',num2str(SR),' Hz'])
else  
plot(t,y,'b')
xlabel(['(a)  An analog ',num2str(Freq),' Hz sine wave.  X-axis = Time, sec'])
end

ylabel(['Volts'])
axis([0 inf 1.1*min(y) 1.1*max(y)])

subplot(312)

hold on

maxy = max(abs(y));

alty = y - min(y);

CodingMethod = CodeMeth;

if CodingMethod==1 % Sign + Mag
   LSBAmt = 1*maxy/((2^(BitsQ-1))-1);
newY = y; 

preOffsetInput = Bias; % vals 1,2 none, val 3 one-half LSB
if preOffsetInput==0
   theOffset = 0;
elseif preOffsetInput==1
   theOffset = LSBAmt/2;
   PosYs = find(newY>=0);
   newY(PosYs) = newY(PosYs)+theOffset;
   NegYs = find(newY<0);
   newY(NegYs) = newY(NegYs)-theOffset;
end

valNormaOrNo = PlotType;  % val 0 plot as multiples of LSB, or 1 to plot in volts

if valNormalOrNo == 0 

reconquantabsy = fix(newY/LSBAmt);
BinOut = LVxADCPosNegQuants(reconquantabsy);
fullamprecon = reconquantabsy; 
fullnoise = fullamprecon - (y/LSBAmt); 
hold on
stem(0:1:length(fullamprecon)-1,fullamprecon,'bo');
plot(0:1:length(fullamprecon)-1,y/LSBAmt,'k')
xlabel(['(b)  Quantized to ',num2str(BitsQ),' bits: 1 sign, ',num2str(BitsQ-1),' mag bit(s).  X-axis = Sample; LSB = ',num2str(LSBAmt),' Volts'])
ylabel(['Multiples/LSB'])
axis([0 length(fullamprecon)-1 1.2*min(y/LSBAmt) 1.2*max(y/LSBAmt)])
else
reconquantabsy = fix(newY/LSBAmt); 
BinOut = LVxADCPosNegQuants(reconquantabsy);
fullamprecon = reconquantabsy*LSBAmt;
fullnoise = fullamprecon - y;  
hold on
stem(0:1:length(fullamprecon)-1,fullamprecon,'bo');
plot(0:1:length(fullamprecon)-1,y,'k')
xlabel(['(b)  Quantized to ',num2str(BitsQ),' bits: 1 sign, ',num2str(BitsQ-1),' mag bit(s).  X-axis = Sample; LSB = ',num2str(LSBAmt),' Volts'])
ylabel(['Volts'])
axis([0 length(fullamprecon)-1 1.1*min(y) 1.1*max(y)])    
end

elseif CodingMethod==2 % Offset
   LSBAmt = 1*max(alty)/(2^(BitsQ)-1);
   
preOffsetInput = Bias; 
if preOffsetInput==0
   theOffset = 0;
elseif preOffsetInput==1
   theOffset = LSBAmt/2;
end

if valNormalOrNo == 0 
reconquantabsy = fix((alty+theOffset)/LSBAmt);
BinOut = LVxADCOffset(reconquantabsy,BitsQ);
fullamprecon = reconquantabsy; 
fullnoise = fullamprecon - (alty/LSBAmt);  
hold on
stem(0:1:length(fullamprecon)-1,fullamprecon,'bo');
plot(0:1:length(fullamprecon)-1,alty/LSBAmt,'k')
xlabel(['(b)  Quantized to ',num2str(BitsQ),' bits, Offset Method.  X-axis = Sample; LSB = ',num2str(LSBAmt),' Volts'])
ylabel(['Multiples/LSB'])
axis([0 length(fullamprecon)-1 1.2*min(alty/LSBAmt) 1.2*max(alty/LSBAmt)])
else
reconquantabsy = fix((alty+theOffset)/LSBAmt);
BinOut = LVxADCOffset(reconquantabsy,BitsQ);
fullamprecon = reconquantabsy*LSBAmt;
fullnoise = fullamprecon - alty;  
hold on
stem(0:1:length(fullamprecon)-1,fullamprecon,'bo');
plot(0:1:length(fullamprecon)-1,alty,'k')
xlabel(['(b)  Quantized to ',num2str(BitsQ),' bits, Offset Method.  X-axis = Sample; LSB = ',num2str(LSBAmt),' Volts'])
ylabel(['Volts'])
axis([0 length(fullamprecon)-1 1.1*min(alty) 1.1*max(alty)])    
end
end

subplot(313)
thenoise = fullnoise;

if valNormalOrNo == 0
stem(0:1:length(thenoise)-1,thenoise,'bo');
xlabel(['(c)  Quantization Noise. X-Axis = Sample; LSB = ',num2str(LSBAmt),' Volts'])
ylabel(['Mults/LSB'])
axis([0 length(thenoise)-1 1.2*min(thenoise) 1.2*max(thenoise)])
else
stem(0:1:length(thenoise)-1,thenoise,'bo');
xlabel(['(c)  Quantization Noise. X-Axis = Sample; LSB = ',num2str(LSBAmt),' Volts'])
ylabel(['Volts'])
axis([0 length(thenoise)-1 1.1*min(thenoise) 1.1*max(thenoise)])    
end

return

function  BinaryOutput = LVxADCPosNegQuants(NumsToDig)
%function  BinaryOutput = LVxADCPosNegQuants(NumsToDig)
%Returns a matrix of binary numbers in sign plus
%magnitude format constituting the binary conversion of the decimal numbers
%input as the vector NumsToDig
% Test call:
% [BinOut] = LVxADCPosNegQuants([17 7 -100  77  -63  98  21  -35])

NumsToDig = real(NumsToDig);
LargestMagnitude = max(abs(NumsToDig));      
n = 1;
while ((2^n)-1)<LargestMagnitude
   n = n+1;
end
NumberOfBits = n + 1; 
[nRows,nCols] = size(NumsToDig);
if nRows==1 & nCols>1
   howManySamps = nCols; 
elseif nRows > 1 & nCols==1
   NumsToDig = NumsToDig';
   howManySamps = nRows;
elseif nRows==1 & nCols==1
   howManySamps=1;
elseif nRows>1 & nCols>1
	Comment = 'Improper format; ending procedure'
   return
end
NegIndices = find(NumsToDig<0);
NumsToDig = abs(NumsToDig);

WtVec(1,1:NumberOfBits-1) = 2.^(0:1:NumberOfBits-2); % vectorized 

OutputLessSignBit = DigitizePosNums(NumsToDig,NumberOfBits-1);
BinaryOutput(1:howManySamps,1) = zeros(howManySamps,1);
BinaryOutput(1:howManySamps,2:NumberOfBits) = OutputLessSignBit;
BinaryOutput(NegIndices,1) = ones(length(NegIndices),1);
return
function [OutputMat] = DigitizePosNums(xNumsToDig,NumBits) 
xWtVec = 2.^(NumBits-1:-1:0);
xhowManySamps = length(xNumsToDig);
OutputMat = zeros(xhowManySamps,NumBits);
for BitCtr = 1:1:NumBits 
    OutputMat(1:xhowManySamps,BitCtr) = 1; 
    TestMat = OutputMat*xWtVec' - xNumsToDig';   
    [i,j] = find(TestMat>0);
    OutputMat(i,BitCtr) = 0;
end
return

function FinOut = LVxADCOffset(NumsToDig,NumBits)

% FinOut = LVxADCOffset([-255,-254,-253,253,254,256],4)
% FinOut = LVxADCOffset([-255,-254,-253,253,254,256],9)
% FinOut = LVxADCOffset([-3,-2,-1,0,1,2,3,4,5],4)

NumsToDig = real(NumsToDig);
theMin = min(NumsToDig);
theMax = max(NumsToDig);
NumsToDig = NumsToDig - theMin; % normalized so lowest value is zero

[nRows,nCols] = size(NumsToDig);
if nRows==1 & nCols>1
   howManySamps = nCols; 
elseif nRows > 1 & nCols==1
   NumsToDig = NumsToDig';
   howManySamps = nRows;
elseif nRows==1 & nCols==1
   howManySamps=1;
elseif nRows>1 & nCols>1
%	Comment = 'Improper format; ending procedure'
   return
end

NumsToDig = abs(NumsToDig);
WtVec(1,1) = 1;
for ctr = 1:1:NumBits-1
    WtVec(1,ctr+1) = 2^ctr;
end
Output = DigitizePosNums2(NumsToDig,NumBits,howManySamps,WtVec);
FinOut = Output;
return

function [OutputMat] = DigitizePosNums2(xNumsToDig, NumBits, xhowManySamps,xWtVec) 
BitWtMat = zeros(xhowManySamps,NumBits);
for BitCtr = NumBits:-1:1 
    BitWtMat(1:xhowManySamps,BitCtr) = 1; 
    TestMat = BitWtMat*xWtVec' - xNumsToDig';   
    [i,j] = find(TestMat>0);
    BitWtMat(i,BitCtr) = 0;
 end
OutputMat = fliplr(BitWtMat);






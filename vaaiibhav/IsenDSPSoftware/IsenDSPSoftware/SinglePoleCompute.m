function SinglePoleCompute
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
global FirstDelayOutput 
global FirstStageFeedbackValue FirstDelayInput
global FirstFeedbackGain RealInputAxis
global  FirstOutputImagAxis
global FirstOutputRealAxis FirstRealPole
global FirstImagPole 
global CompModeSelect100 SeqLengthSet TestSignalType
global UCPolePlot 

ftsz1 = 13;
ftsz2 = 13;
ftsz3 = 13;
RTP = str2num(get(FirstRealPole,'String'));
ITP = str2num(get(FirstImagPole,'String'));

TheAngInRadians = (2*pi*(ITP/360));

stemtoplotlim = 66;

SigDig = 2;

ThePole = RTP*( cos(TheAngInRadians) + j*sin(TheAngInRadians) );
if abs(real(ThePole))<10^-4
ThePole = RTP*(j*sin(TheAngInRadians) );   
elseif abs(imag(ThePole))<10^-4
ThePole = RTP*(cos(TheAngInRadians));  
end

set(UCPolePlot,'xdata',real(ThePole));
set(UCPolePlot,'ydata',imag(ThePole));

set(FirstFeedbackGain,'string',num2str(ThePole,SigDig))
%set(SecondFeedbackGain,'String',num2str(CCThePole));
%set(FirstFeedbackGain,'String',num2str(ThePole));

GetSeqLen = get(SeqLengthSet,'Value');
if (GetSeqLen>1)
   if GetSeqLen==2
        SR = 4;
    elseif GetSeqLen==3
        SR = 8;
	elseif GetSeqLen==4
        SR = 16;  
	elseif GetSeqLen==5
        SR = 32;  
	elseif GetSeqLen==6
        SR = 64;
    elseif GetSeqLen==7
        SR = 128;
 	elseif GetSeqLen==8
        SR = 256;
    elseif GetSeqLen==9
        SR = 512;
    elseif GetSeqLen==10
        SR = 1024;     
   end
else
   SR = 32;
end

t = 0:1/SR:1-1/SR;

TST = get(TestSignalType ,'Value');
if (TST>1)
   if TST==2
   xInput(1,1:length(t))=zeros; %use this to form an impulse 
   xInput(1,1) = 1;
     elseif TST==3
        xInput(1,1:length(t)) = ones;
        %xInput(1,1:1) = 0;
     elseif TST==4
      xInput(1,1:length(t)) = chirp(t,1,1,SR/2);
	  end
else
   xInput(1,1:length(t)) = chirp(t,1,1-1/SR,SR/2);
end

b = 1;
a = [1,-(ThePole)];
FirstAnswer = filter(b,a,xInput);

RealMax = max((max(abs((FirstAnswer)))));
ImagMax = RealMax;

Output1(1,1:SR) = zeros;
Output1(1) = xInput(1);

CMS = get(CompModeSelect100 ,'Value');
if (CMS>1)
   if CMS==2
      CompMode = 1;  
		elseif CMS==3
			CompMode = 2; 
   	elseif CMS==4
			CompMode = 3; 
	end
else
   CompMode = 3;
end

if CompMode == 3
  subplot(RealInputAxis)
if SR< stemtoplotlim
   hndPlot = stem(xInput(1:SR),'ko');
   set(hndPlot,'markersize',3)
else
   plot(xInput(1:SR),'b')
end
ylabel('Real','fontsize',ftsz2)
xlabel('(a) Sample Number of Input','fontsize',ftsz1)
set(gca,'fontsize',ftsz3)
axis([0 SR+1 -1.5 1.5])

subplot(FirstOutputRealAxis)
if SR< stemtoplotlim
   hndPlot = stem(real(FirstAnswer(1:SR)),'ko');
   set(hndPlot,'markersize',3)
else
   plot(real(FirstAnswer(1:SR)),'k')   
end
ylabel('Real','fontsize',ftsz2)
xlabel('(b) Sample Number of Real Output','fontsize',ftsz1)
set(gca,'fontsize',ftsz3)
axis( [ 0  SR+1  -1.2*(RealMax)  1.2*(RealMax) ])

subplot(FirstOutputImagAxis)
if SR< stemtoplotlim
   hndPlot = stem(imag(FirstAnswer(1:SR)),'ko');
   set(hndPlot,'markersize',3)
else
   plot(imag(FirstAnswer(1:SR)),'k')   
end
ylabel('Imag','fontsize',ftsz2)
xlabel('(c) Sample Number of Imaginary Output','fontsize',ftsz1)
set(gca,'fontsize',ftsz3)
axis( [ 0  SR+1  -1.2*(ImagMax)  1.2*(ImagMax) ])

return
end

Output1(1) = xInput(1);

for CurDataPtr1 = 2:1:SR
   FirstFB = ThePole*Output1(CurDataPtr1 - 1);
   %FirstFB = ZeroSmallValues(FirstFB);
   
   stFB1 = num2str(FirstFB,SigDig);
   LenstFB1 = length(stFB1);
   stFB1 = [blanks(13-LenstFB1),stFB1];
   
   set(FirstStageFeedbackValue,'String',stFB1)
   
   Output1(CurDataPtr1) = xInput(CurDataPtr1) + FirstFB;
   
   xOutput1 = ZeroSmallValues(Output1(CurDataPtr1));
   
   stDelIn1 = num2str(xOutput1,SigDig);
   LenstDelIn1 = length(stDelIn1);
   stDelIn1 = [blanks(13-LenstDelIn1),stDelIn1];
   
   xxxOutput1 = ZeroSmallValues(Output1(CurDataPtr1-1));
   stDelOut1 = num2str(xxxOutput1,SigDig);
   LenstDelOut1 = length(stDelOut1);
   stDelOut1 = [blanks(13-LenstDelOut1),stDelOut1];
   
   set(FirstDelayInput,'String',stDelIn1)
   set(FirstDelayOutput,'String',stDelOut1)
     
subplot(RealInputAxis)
if SR< stemtoplotlim
   hndPlot = stem(xInput(1:SR),'ko');
   set(hndPlot,'markersize',3)
else
   plot(xInput(1:SR),'k')
end
ylabel('Real','fontsize',ftsz2)
xlabel('(a) Sample Number of Input','fontsize',ftsz1)
set(gca,'fontsize',ftsz3)
axis([0 SR+1 -1.5 1.5])

subplot(FirstOutputRealAxis)
if SR< stemtoplotlim
   hndPlot = stem(real(Output1(1:SR)),'ko');
   set(hndPlot,'markersize',3)
else
   plot(real(Output1(1:SR)),'k')   
end
ylabel('Real','fontsize',ftsz2)
xlabel('(b) Sample Number of Real Output','fontsize',ftsz1)
set(gca,'fontsize',ftsz3)
axis( [ 0  SR+1  -1.2*(RealMax)  1.2*(RealMax) ])

subplot(FirstOutputImagAxis)
if SR< stemtoplotlim
   hndPlot = stem(imag(Output1(1:SR)),'ko');
   set(hndPlot,'markersize',3)
else
   plot(imag(Output1(1:SR)),'k')   
end
ylabel('Imag','fontsize',ftsz2)
xlabel('(c) Sample Number of Imaginary Output','fontsize',ftsz1)
set(gca,'fontsize',ftsz3)
axis( [ 0  SR+1  -1.2*(ImagMax)  1.2*(ImagMax) ])

if CompMode == 1 
      pause(0.1)
elseif CompMode == 2    
      pause(1.25)
end
   
end

function FixedNumber = ZeroSmallValues(NumberArg)

FixedNumber = NumberArg;

if abs(real(NumberArg))<0.0001
   FixedNumber = 0 + imag(NumberArg);
end

if abs(imag(NumberArg))<0.0001
   FixedNumber = real(FixedNumber) + 0*i;
end




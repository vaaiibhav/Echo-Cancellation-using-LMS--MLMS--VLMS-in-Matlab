function SuccApproxCompute
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
global hndMSBProdValBox hndSecMSBProdValBox hndThirdMSBProdValBox
global hndFourthMSBProdValBox selStepTypeSAR hndInputSigAxis
global hndTwoToThreeSAR hndTwoToTwoSAR hndTwoToOneSAR hndTwoToZeroSAR
global hndLSBAmpSet hndSetPeakValTestSig hndCommentBox hndFigSuccApprox
global hndOutputDisplayAxis hndSumOut hndSRSet hndWFSel
global hndInputSigValBox hndDACOutToCompareIn hndConversionBox
global MovingArrowLeadLine MovingArrowHead hndBiasSet
global TwoToThreeDACSource TwoToTwoDACSource TwoToOneDACSource TwoToZeroDACSource

figure(hndFigSuccApprox)
subplot(hndOutputDisplayAxis)
cla

ftsz = [14];

set(hndTwoToThreeSAR,'String','0')
set(hndTwoToTwoSAR,'String','0')
set(hndTwoToOneSAR,'String','0')
set(hndTwoToZeroSAR,'String','0')

set(hndMSBProdValBox,'String','0')
set(hndSecMSBProdValBox,'String','0')
set(hndThirdMSBProdValBox,'String','0')
set(hndFourthMSBProdValBox,'String','0')

set(hndSumOut,'String','0')

POutputBias = get(hndBiasSet,'Value');
if POutputBias==1|POutputBias==2
   OutputBias = 2;
elseif POutputBias==3
    OutputBias=3;
else
   OutputBias = 2;
end

PSR = get(hndSRSet,'Value');
if PSR>1
   if PSR==2
      SR = 8;
   elseif PSR==3
      SR=16;
   elseif PSR==4
      SR=32;
   elseif PSR==5
      SR=64;
   elseif PSR==6
      SR=128;      
   elseif PSR==7
      SR=256;
   end
else
   SR = 32;
end

AnsSet(1,SR) = zeros;

PCompMode = get(selStepTypeSAR,'Value');
if PCompMode>1
   CompMode=PCompMode;
else
   CompMode=2;
end

if CompMode==4
   SR = 8;
end

pkInputAmp = get(hndSetPeakValTestSig,'Value');

if pkInputAmp > 1
   if pkInputAmp==2
      InputAmp = 1;
   elseif pkInputAmp==3
      InputAmp = 2;
   elseif pkInputAmp==4
      InputAmp = 3;
   elseif pkInputAmp==5
      InputAmp = 7;
   elseif pkInputAmp==6
      InputAmp = 15;
   elseif pkInputAmp==7
      InputAmp = 31;
   elseif pkInputAmp==8
      InputAmp = 63;
   end
   
else
   InputAmp = 15;
end

%InputAmp = InputAmp/2;

amtLSB = get(hndLSBAmpSet,'Value');
if amtLSB>1
   if amtLSB==2
      LSBAmp = 1;
    elseif amtLSB==3
      LSBAmp = 2;
   elseif amtLSB==4
      LSBAmp = 3;
   elseif amtLSB==5
      LSBAmp = 4;   
   elseif amtLSB==6
      LSBAmp = 5;        
   else
      LSBAmp = 1;
   end
else
   LSBAmp = 1;

end

set(TwoToThreeDACSource,'String',num2str(8*LSBAmp))
set(TwoToTwoDACSource,'String',num2str(4*LSBAmp))
set(TwoToOneDACSource,'String',num2str(2*LSBAmp))
set(TwoToZeroDACSource,'String',num2str(LSBAmp))

%SR = 64;
t = 0:1/SR:1-1/SR;
SelWF = get(hndWFSel,'Value');
if SelWF>1
   if SelWF==2
      TestSig = InputAmp*sin(2*pi*t);
   elseif SelWF==3
      TestSig = InputAmp*cos(2*pi*t);
   elseif SelWF==4
      TestSig = InputAmp*sin(2*pi*t*2);
   elseif SelWF==5
      TestSig = InputAmp*cos(2*pi*t*2);
   elseif SelWF==6        
		TestSig = InputAmp*sin(2*pi*t) + 0.5*InputAmp*cos(2*pi*t*2.5);
	elseif SelWF==7        
		TestSig = InputAmp*sin(2*pi*t) + 0.5*InputAmp*cos(2*pi*t*3.5);
   end
else
      TestSig = InputAmp*sin(2*pi*t) + 0.5*InputAmp*cos(2*pi*t*2.5);
end
   
TestSig = TestSig -  min(TestSig);
TestSig = TestSig*(InputAmp/max(TestSig));
OrigTestSig = TestSig;

if OutputBias==2
elseif OutputBias==3
    TestSig = TestSig + LSBAmp/2;  
end

subplot(hndInputSigAxis)

plot(0:1:length(t)-1,OrigTestSig)
axis([0 length(t)-1  0 1.2*InputAmp])

xlabel(['(a) Input Sample'],'fontsize',ftsz)
ylabel(['Volts'],'fontsize',ftsz)

subplot(hndOutputDisplayAxis)

plot(0:1:length(t)-1,OrigTestSig,'k:')
%axis([0 length(t)-1  -1.2*InputAmp 1.2*InputAmp])
axis([0 length(t)-1  0 1.2*InputAmp])
xlabel(['(b) Input Signal (Dashed); Quantized Samples (Stem Plot); X-Axis = Sample No.'],'fontsize',ftsz)
ylabel(['Volts'],'fontsize',ftsz)

finalConv = [num2str(0),' ',num2str(0),' ',num2str(0),' ',num2str(0)];   

for GrandCtr = 0:1:length(TestSig)-1

   subplot(hndInputSigAxis)
   	cla
		hold on
      plot(0:1:length(t)-1,TestSig)
    %  plot(GrandCtr,TestSig(GrandCtr+1),'ro')
      stem(GrandCtr,TestSig(GrandCtr+1),'ro')
		axis([0 length(t)-1  0 1.2*InputAmp])
        set(gca,'fontsize',18)
      xlabel(['(a) Input Signal'],'fontsize',ftsz)
      ylabel(['Volts'],'fontsize',ftsz)
      set(hndInputSigValBox,'String',num2str(TestSig(GrandCtr+1),3)  )
           
      curSigIn = TestSig(GrandCtr+1);
      
curvalTwoToThree = 0;
curvalTwoToTwo = 0;
curvalTwoToOne = 0;
curvalTwoToZero = 0;

MSBProd = 0;
SecMSBProd = 0;
ThirdMSBProd = 0;
FourthMSBProd = 0;

set(hndTwoToThreeSAR,'String','0')
set(hndTwoToTwoSAR,'String','0')
set(hndTwoToOneSAR,'String','0')
set(hndTwoToZeroSAR,'String','0')

set(hndMSBProdValBox,'String','0')
set(hndSecMSBProdValBox,'String','0')
set(hndThirdMSBProdValBox,'String','0')
set(hndFourthMSBProdValBox,'String','0')

for QuantCtr = 3:-1:0
   
XOffset = (3 - QuantCtr)*12;
   
set(MovingArrowLeadLine,...
'XData',[-31.25000000000001 -27.71739130434783 -27.71739130434783 -18.47826086956522+XOffset -18.47826086956522+XOffset], ...
'YData',[20.73170731707315 20.73170731707315 89.43089430894307 89.43089430894307 80.08130081300811]);
 
set(MovingArrowHead,...
'Vertices',[-19.2783+XOffset 84.0813; -18.4783+XOffset  80.0813; -17.6783+XOffset  84.0813; -19.2783+XOffset 84.0813]);

set(hndCommentBox,'String',['    Set the 2^{',num2str(QuantCtr,1),'} bit to 1'],'fontsize',ftsz)

if QuantCtr==3
set(hndTwoToThreeSAR,'String','1')
curvalTwoToThree = 1;     
set(hndMSBProdValBox,'String','8')   
elseif QuantCtr==2
set(hndTwoToTwoSAR,'String','1')
curvalTwoToTwo = 1;     
set(hndSecMSBProdValBox,'String','4')    
elseif QuantCtr==1
set(hndTwoToOneSAR,'String','1')
curvalTwoToOne = 1;     
set(hndThirdMSBProdValBox,'String','2')    
elseif QuantCtr==0
set(hndTwoToZeroSAR,'String','1')
curvalTwoToZero = 1;     
set(hndFourthMSBProdValBox,'String','1')    
end
  
  MSBProd = curvalTwoToThree*8*LSBAmp;
  set(hndMSBProdValBox,'String',num2str(MSBProd))
  
  SecMSBProd = curvalTwoToTwo*4*LSBAmp;
  set(hndSecMSBProdValBox,'String',num2str(SecMSBProd))
  
  ThirdMSBProd = curvalTwoToOne*2*LSBAmp;
  set(hndThirdMSBProdValBox,'String',num2str(ThirdMSBProd))
  
  FourthMSBProd = curvalTwoToZero*1*LSBAmp;
  set(hndFourthMSBProdValBox,'String',num2str(FourthMSBProd))
  
  curDACOutput = sum(MSBProd + SecMSBProd + ThirdMSBProd + FourthMSBProd);
  
set(hndSumOut,'String',num2str(curDACOutput,2))
set(hndDACOutToCompareIn,'String',num2str(curDACOutput,2))

%LogicAns = CompareVals(curDACOutput,curSigIn);
   if curSigIn-curDACOutput < 0
      valCompare=0;
   else
      valCompare=1;
   end
   
if valCompare==0
      
    		if CompMode==2
   		pause(0.01) 
			elseif CompMode==3
   		pause(0.7)
			elseif CompMode==4
   		pause
			else
   		pause(0.01)
			end
       
   if QuantCtr==3
      curvalTwoToThree=0;
      set(hndTwoToThreeSAR,'String','0')
      set(hndCommentBox,'String',['Reset the 2^{',num2str(QuantCtr),'} bit to 0'],'fontsize',ftsz)
   elseif QuantCtr==2
      curvalTwoToTwo=0;
      set(hndTwoToTwoSAR,'String','0')
      set(hndCommentBox,'String',['Reset the 2^{',num2str(QuantCtr),'} bit to 0'],'fontsize',ftsz)
   elseif QuantCtr==1
   curvalTwoToOne=0;
   set(hndTwoToOneSAR,'String','0')
   set(hndCommentBox,'String',['Reset the 2^{',num2str(QuantCtr),'} bit to 0'],'fontsize',ftsz)
   elseif QuantCtr==0
   curvalTwoToZero=0;
   set(hndTwoToZeroSAR,'String','0')
   FourthMSBProd = 0; %curvalTwoToZero*1*LSBAmp;
  	set(hndFourthMSBProdValBox,'String','0'); %num2str(FourthMSBProd)) 
	curDACOutput = sum(MSBProd + SecMSBProd + ThirdMSBProd + FourthMSBProd);
   set(hndSumOut,'String',num2str(curDACOutput,2))
   set(hndDACOutToCompareIn,'String',num2str(curDACOutput,2))
   set(hndTwoToZeroSAR,'String','0')

   set(hndCommentBox,'String',['Reset the 2^{',num2str(QuantCtr),'} bit to 0'],'fontsize',ftsz)
   finalConv = [num2str(curvalTwoToThree),' ',num2str(curvalTwoToTwo),' ',num2str(curvalTwoToOne),' ',num2str(curvalTwoToZero)];   
	end
else % compare was not 0
   if QuantCtr==0
		finalConv = [num2str(curvalTwoToThree),' ',num2str(curvalTwoToTwo),' ',num2str(curvalTwoToOne),' ',num2str(curvalTwoToZero)];
		set(hndConversionBox,'String',['Conversion Value for Sample ',num2str(GrandCtr),' is ',finalConv],'fontsize',ftsz);     
    end  
end
set(hndConversionBox,'String',['Conversion Value for Sample ',num2str(GrandCtr),' is ',finalConv],'fontsize',ftsz);     

if CompMode==2
   pause(0.01) 
elseif CompMode==3
   pause(1)
elseif CompMode==4
   pause
else
   pause(0.01)
end

end % for quantizing loop

subplot(hndOutputDisplayAxis)
%cla
hold on

%if OutputBias==2
   quantans = curDACOutput;
   %else
%   quantans = curDACOutput + LSBAmp/2;
%end

AnsSet(GrandCtr+1) = quantans;

hndPlot = stem(GrandCtr,quantans,'bo');
set(hndPlot,'markersize',4)

axis([0 length(t)-1  0 1.2*InputAmp])
set(gca,'fontsize',18)
%xlabel(['The Input Signal (Dashed); Quantized Samples (Stem Plot); X-Axis = Sample No.'],'fontsize',ftsz)
xlabel(['(b) Input Signal (Dashed); Quantized Samples (Stem Plot); X-Axis = Sample No.'],'fontsize',ftsz)
ylabel(['Volts'],'fontsize',ftsz)


if CompMode==2
   pause(0.01) 
elseif CompMode==3
   pause(1)
elseif CompMode==4
   pause
else
   pause(0.01)
end

end % for sample ctr

%break
scnsize = get(0,'ScreenSize');
pos1 = [0.01*scnsize(3), 0.05*scnsize(4), 0.98*scnsize(3), 0.84*scnsize(4)];

figure(7003)
clf
set(7003,'color',[1 1 1],'position',pos1,'numbertitle','off','name','Results of Quantization')

hold on
stairs(0:1:(length(AnsSet)-1),AnsSet,'k')
plot(0:1:(length(AnsSet)-1),OrigTestSig,'ko')
plot(0:1:(length(AnsSet)-1),OrigTestSig,'k:')
hilim = max(abs(AnsSet));
% xlabel(['Analog Signal (Dashed); Ideal Samples (Circles); Quantized Samps at DAC output (Solid)'],'fontsize',[12])
set(gca,'fontsize',18)
xlabel(['Sample'],'fontsize',[18])
ylabel(['Amplitude, Decimal Numbers'],'fontsize',[18])
axis([0 SR-1 -1 1.2*hilim])

hold off

figure(hndFigSuccApprox)

%figure(7004)
%clf
%set(7004,'color',[1 1 1])

%semilogy(abs(fft(TestSig,1024)),'k')
%hold on

%semilogy(abs(fft(AnsSet,1024)),'k:');

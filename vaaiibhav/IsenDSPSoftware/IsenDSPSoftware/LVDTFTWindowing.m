function LVDTFTWindowing(SR)
%function LVDTFTWindowing(SR)
%
%SR is the length in samples of each of several windows to be displayed;
%The window, the window imposed on a sequence at the half-band
%frequency, and the spectrum of the latter are all displayed.
%
%A typical call might be 
% LVDTFTWindowing(32)

% Copyright 2007 by F. W. Isen

dispszXL = 15;

figure(7766)

t = 0:1/SR:1-1/SR;
countlim = 15;

for ctr = 1:1:countlim

  if ctr==1
      kfac = '';
      theWin = 'rectwin';
      stTheWin = 'Rectangular';
    elseif ctr==2
           kfac = '';
           theWin = 'hann';
           stTheWin = 'Hanning';
    elseif ctr==3
           kfac = '';
           theWin = 'hamming';
           stTheWin = 'Hamming';
    elseif ctr==4
           kfac = '';
           theWin = 'triang';
           stTheWin = 'Triangle';
    elseif ctr==5
           kfac = '';
           theWin = 'blackman';
           stTheWin = 'Blackman';
    elseif ctr==6
           kfac = '';
           theWin = 'bartlett';
           stTheWin = 'Bartlett';
    elseif ctr==7
       kfac = '2';
       theWin = 'kaiser';
       stTheWin = 'Kaiser';
    elseif ctr==8
       kfac = '3';
       theWin = 'kaiser';
       stTheWin = 'Kaiser';
    elseif ctr==9
       kfac = '4';
       theWin = 'kaiser';
       stTheWin = 'Kaiser';
    elseif ctr==10
       kfac = '5';
       theWin = 'kaiser';
       stTheWin = 'Kaiser';
    elseif ctr==11
       kfac = '20';
       theWin = 'chebwin';
       stTheWin = 'Chebwin';
    elseif ctr==12
       kfac = '30';
       theWin = 'chebwin';
       stTheWin = 'Chebwin';
    elseif ctr==13
       kfac = '40';
       theWin = 'chebwin';
       stTheWin = 'Chebwin';
    elseif ctr==14
       kfac = '50';
       theWin = 'chebwin';
       stTheWin = 'Chebwin';
    elseif ctr==15
       kfac = '60';
       theWin = 'chebwin';
       stTheWin = 'Chebwin';
    end 
   
if ctr==1
    TDSeq1 = rectwin(SR);
elseif ctr==2
    TDSeq1 = hanning(SR);
elseif ctr==3
    TDSeq1 = hamming(SR);
elseif ctr==4
    TDSeq1 = triang(SR);
elseif ctr==5
    TDSeq1 = blackman(SR);
elseif ctr==6
     TDSeq1 = bartlett(SR);               
elseif ctr==7
   TDSeq1 = kaiser(SR,str2num(kfac));
elseif ctr==8
   TDSeq1 = kaiser(SR,str2num(kfac));
elseif ctr==9
   TDSeq1 = kaiser(SR,str2num(kfac));
elseif ctr==10
   TDSeq1 = kaiser(SR,str2num(kfac));
elseif ctr>10
   TDSeq1 = chebwin(SR,str2num(kfac));
end

TDSeq1 = TDSeq1';
Freq = SR/4;
TDSeq2 = sin(2*pi*t*Freq); 
plotlim1 = 1.3*max(abs(TDSeq1));

subplot(221)

stem(0:1:length(TDSeq1)-1,TDSeq1);
ylabel(['Amp'])
if ctr==countlim
xlabel(['(a)  Window: ',stTheWin,' ',kfac])
else   
xlabel(['(a)  Window: ',stTheWin,' ',kfac])
end
axis([0, length(TDSeq1)-1, 0,  plotlim1])

plotlim2 = 1.75*max(abs(TDSeq2));

subplot(222)
stem(0:1:length(TDSeq2)-1,TDSeq2);
ylabel(['Amp'])
xlabel(['(b) Sequence'])
axis([0, length(TDSeq2)-1,  -plotlim2,  plotlim2])

ProdTDSeq1TDSeq2 = TDSeq1.*TDSeq2;
plotlim3 = 1.75*max(abs(ProdTDSeq1TDSeq2));

subplot(223)
stem(0:1:length(TDSeq2)-1,ProdTDSeq1TDSeq2);
ylabel(['Amp'])
xlabel(['(c) (Window.*Sequence) (The Windowed Sequence)'])
axis([0, length(TDSeq2)-1,  -plotlim3,  plotlim3])

dtftProd = abs(fft(ProdTDSeq1TDSeq2,32*SR));
dtftProd = dtftProd/max(abs(dtftProd));

plotlim4 = 1.15*max(dtftProd);

subplot(224)

NyqL = length(dtftProd)/2;

plot(0:1/NyqL:(NyqL-1)/NyqL,20*log10(dtftProd(1,1:length(dtftProd)/2)+eps));
grid on
ylabel(['Mag'])
xlabel(['(d)  Normalized Spectrum of Windowed Sequence'])
axis([0, (NyqL-1)/NyqL,  -80,  10])

pause; %(1)
end




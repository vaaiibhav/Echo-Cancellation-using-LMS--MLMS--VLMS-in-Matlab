function LVDTFTWindowsOnly(SR)
% function LVDTFTWindowsOnly(SR)
%
% Displays a series of windows and their DTFTs
% SR is the length of each window to be displayed.
%
% A typical call is:
%
% LVDTFTWindowsOnly(32)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

dispszXL = 12;

figure(7767)
clf

countlim = 15;

for ctr = 1:1:countlim
   
   if ctr==1
      kfac = '';
      theWin = 'rectwin';
      stTheWin = 'Rectwin';
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
    TDSeq1 = hann(SR);
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

Freq = 6*(SR/32);

plotlim1 = 1.3*max(abs(TDSeq1));

subplot(211)
stem(0:1:length(TDSeq1)-1,TDSeq1,'bo');
ylabel(['Amplitude'])
if ctr==countlim
xlabel(['(a)  The Window: ',stTheWin,' ',kfac,':  (Last Window)'])
else   
xlabel(['(a)  The Window: ',stTheWin,' ',kfac,':  (Press Any Key to See the Next Window)'])
end
axis([0 length(TDSeq1)-1 0  plotlim1])

dtftProd = (abs(fft(TDSeq1,16*SR)))+ eps;
dtftProd = 20*log10(dtftProd/(max(dtftProd)));

subplot(212)

NyqL = length(dtftProd)/2;
plot(0:1/NyqL:(NyqL-1)/NyqL,dtftProd(1,1:length(dtftProd)/2),'b');
ylabel(['Mag, dB'])
xlabel(['(b)  Normalized Spectrum of the ',stTheWin,' Window Above'])
axis([0 (NyqL-1)/NyqL -120 5])

pause

end

function LVCircularConvolution(TDSeq1,TDSeq2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
figure(10987)

Signal = TDSeq1;
TDHilbert = TDSeq2;

subplot(321)
stem(real(Signal),'bo');
xlabel(['(a) Sample, 1st Seq'])
ylabel(['Real'])
axis([0 (length(Signal)+1) 1.2*min(real(Signal))  1.2*max(real(Signal))])

subplot(322)
stem(imag(Signal),'bo');
xlabel(['(b) Sample, 1st Seq'])
ylabel(['Imag'])
axis([0 (length(Signal)+1) 1.2*min(real(Signal))  1.2*max(real(Signal)) ])

subplot(323)
stem(real(TDHilbert),'bo');
xlabel(['(c) Sample, 2nd Seq'])
ylabel(['Real'])
axis([0 (length(TDHilbert)+1) 1.2*1.2*min(min([real(TDHilbert) imag(TDHilbert)]))   1.2*max(max([real(TDHilbert) imag(TDHilbert)])) ])


subplot(324)
stem(imag(TDHilbert),'bo');
xlabel(['(d) Sample, 2nd Seq'])
ylabel(['Imag'])
axis([0 (length(TDHilbert)+1) 1.2*min([min(real(TDHilbert)) min(imag(TDHilbert)) ] )   1.2*max([max(real(TDHilbert))  max(imag(TDHilbert))])    ])
% do circ con of Signal and TDHilbert (complex) and plot real and imag
% parts on lower axes

% CirCon = ifft(fft(Signal).*fft(TDHilbert));

N = length(TDHilbert);
n= 1:1:N;
for m = 1:1:N
CirCon(m) = sum(Signal(n).*(TDHilbert(mod(m - n, N) +1) ) );
end

subplot(325)
stem(real(CirCon),'bo');
xlabel(['(e) Sample, Circ Conv'])
ylabel(['Real'])
axis([0 (length(CirCon)+1) 1.2*min(min([real(CirCon) imag(CirCon)]))   1.2*max(max([real(CirCon) imag(CirCon)])) ])

subplot(326)
stem(imag(CirCon),'bo');
xlabel(['(f) Sample, Circ Conv'])
ylabel(['Imag'])
axis([0 (length(CirCon)+1) 1.2*min( [-1  min(imag(CirCon)) ])   1.2*max([ 1 max(imag(CirCon)) ] )  ])
function LVDemoRemez(wp,ws,As,Rp,L,PassbandType,NormXFr,EqualWt)
% function LVDemoRemez(wp,ws,As,Rp,L,PassbandType,NormXFr,EqualWt)
% wp, ws, and NormXFr are all normalized frequencies, lying between 0 and
% 1.0
% wp is the passband edge for a lowpass or highpass filter
% ws is the stopband edge for a lowpass or highpass filter
% As is the desired stopband attenuation in dB; if using equal weighting, 
% pass as the empty matrix []
% Rp is the desired passband ripple in dB; if using equal weighting, 
% pass as the empty matrix []
% NormXFr is a vector of normalized extremal frequencies
% PassbandType = 1 gives lowpass, PassbandType=2 gives highpass
% EqualWt = 1 gives equal weight to passband and stopband ripple
%
% Sample calls:
%
% LVDemoRemez(0.45,0.55,40,0.5,33,1,[],1) % u
% LVDemoRemez(0.45,0.55,40,0.5,9,1,[0,0.2875,0.45,0.55,0.7125,1],1) % u
% LVDemoRemez(0.45,0.55,40,0.5,19,1,[0,0.142,0.27,0.385,0.45,0.55,0.59,0.68,0.78,0.89,1],0) % con
% LVDemoRemez(0.45,0.55,30,0.5,9,1,[0,0.29,0.45,0.55,0.72,1],0) % con
% LVDemoRemez(0.65,0.75,60,0.5,19,1,[0.07,0.21,0.332,0.46,0.57,0.65,0.75,0.
% 775,0.835,0.92,1],0) %con
% LVDemoRemez(0.65,0.75,40,0.5,19,1,[0,0.115,0.234,0.355,0.48,0.585,0.65,0.
% 75,0.795,0.89,1],0) % con
% LVDemoRemez(0.65,0.75,40,0.5,19,1,[0,0.135,0.265,0.385,0.5,0.6,0.65,0.75,
% 0.8,0.895,1],1) % u
% LVDemoRemez(0.65,0.75,50,0.5,19,1,[0,0.115,0.234,0.355,0.47,0.585,0.65,...
% 0.75,0.782,0.873,1],0) % con
% LVDemoRemez(0.65,0.75,55,0.5,19,1,[0,0.115,0.234,0.347,0.466,0.578,0.65,
% 0.75,0.782,0.855,1],0) % con
% LVDemoRemez(0.65,0.75,58,0.5,19,1,[0,0.115,0.234,0.347,0.466,0.578,0.65,0
% .75,0.778,0.843,0.93],0) % con
% LVDemoRemez(0.65,0.75,60,0.5,19,1,[0,0.127,0.28,0.39,0.51,0.578,0.75,0.77
% 4,0.835,0.916,1],0) % con
% LVDemoRemez(0.45,0.55,30,0.5,10,1,[0,0.31,0.45,0.55,0.675,0.888],0) % con
%
% Type II--even length
% LVDemoRemez(0.65,0.8,60,0.5,20,1,[0,0.12,0.23,0.357,0.48,0.585,0.65,0.8,0.825,0.885,0.965],0)
% LVDemoRemez(0.65,0.75,40,0.5,20,1,[0,0.106,0.208,0.315,0.41,0.519,0.6062,
% 0.65,0.75,0.8,0.9063],1) % u
%
% HPF Filters, Type I
% LVDemoRemez(0.75,0.65,60,0.5,19,2,[0,0.125,0.27,0.385,0.5,0.6,0.65,0.75,
% 0.795,0.892,1],1) % u
% LVDemoRemez(0.75,0.65,40,0.5,19,2,[0,0.1,0.22,0.35,0.46,0.56,0.61,0.75, 
% 0.79,0.88,1],0) % u
% LVDemoRemez(0.75,0.65,40,0.5,19,2,[0,0.138,0.269,0.3875,0.5,0.6,0.65,0.75
% ,0.8,0.894,1],1) % u
% LVDemoRemez(0.75,0.65,40,0.5,19,2,[0,0.1063,0.2125,0.3187,0.425,0.525,0.6125,0.65,0.75,0.8125,0.9375],0)
% LVDemoRemez(0.75,0.65,40,0.5,19,2,[0,0.106,0.215,0.3175,0.4237,0.523,0.61,0.65,0.75,0.813,0.945],0)
% LVDemoRemez(0.75,0.65,60,0.5,19,2,[0,0.11,0.225,0.335,0.44,0.54,0.615,0.65,0.75,0.85,1],0)
% LVDemoRemez(0.7,0.6,60,0.5,9,2,[0,0.328,0.52,0.6,0.7,1],0)
% LVDemoRemez(0.7,0.6,60,0.5,9,2,[0,0.24,0.47,0.6,0.7,0.85],1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
LenGrid = 2048;

wc = (wp+ws)/2;

if rem(L,2)==0
    kLim = L/2 -1;
    FiltType = 2;
else
    kLim = (L-1)/2;
    FiltType = 1;
end
if (isempty(Rp)|isempty(As))
    EqualWt=1;
else
    EqualWt = 0;
Rfac= 10^(-Rp/20); 
DeltaP = (1-Rfac)/(1+Rfac);
DeltaS = (1+DeltaP)*10^(-As/20);
end

GridpInd = 0:1:LenGrid;
indWc = 0:1:fix(LenGrid*wc);

if PassbandType==1  % lowpass
    BandLims = [0,wp,ws,1];
    Hdr(1,indWc+1) = 1;
    Hdr(1,max(indWc)+1:LenGrid+1) = 0;
    W(1,1:length(indWc)) = DeltaS/DeltaP;
    W(1,length(indWc)+1:LenGrid+1) = 1;
elseif PassbandType==2 % highpass
    BandLims = [0,ws,wp,1];
    Hdr(1,indWc+1) = 0;
    Hdr(1,max(indWc)+1:LenGrid+1) = 1;
    W(1,1:length(indWc)) = 1;
    W(1,length(indWc )+1:LenGrid+1) = DeltaS/DeltaP;
else
    error('PassbandType must be passed as 1 for lowpass or 2 for highpass')
end

if PassbandType==2 & FiltType==2
    error('Type II filter (even length) may not be a highpass filter')
end

if EqualWt==1
    W(1,1:LenGrid+1) = ones;
end

if PassbandType==1 % lowpass
    leftXitionBndLim = wp;
    rightXitionBndlim = ws;
else % highpass
    leftXitionBndLim = ws;
    rightXitionBndlim = wp;
end

% initial guess XFr
NonlinCompres = 0;
% initial guess XFr
if isempty(NormXFr)
NormXFr = [0:1:kLim+1]/(kLim+1);
end

currExtrFreqs = pi*NormXFr;
currNormExtrFreqs = NormXFr

NormFrGrid = [0:1:LenGrid]/LenGrid;
FrGrid = pi*NormFrGrid;

if rem(L,2)==0 % even length, Type II
Q = cos(FrGrid/2);
else
Q = 1; % Type I filter    
end

W = W.*Q;
szW = size(W);
Hdr = Hdr./Q;

if NormXFr(length(NormXFr)) > 1.0
    NormXFr(length(NormXFr)) = 1.0;
end
    
XFrIndonFineGrid = round([NormXFr]*LenGrid) + 1;

if XFrIndonFineGrid(length(XFrIndonFineGrid)) > LenGrid + 1
    XFrIndonFineGrid = XFrIndonFineGrid - (length(XFrIndonFineGrid) - (LenGrid+1));
end

if XFrIndonFineGrid(1) <= 0
    XFrIndonFineGrid(1) = 1;
end

HdrVec = Hdr(1,XFrIndonFineGrid);  % provide this below after finding new extremal freqs
% initial guess comes from code immediately before Grand loop starts
indvalW4XFr = XFrIndonFineGrid; % round(NormXFr*LenGrid)+1;

valW4XFr = W(indvalW4XFr);

WMat(1:kLim+2,1) = 1;

% for k = 0:1:kLim+1
k = 0:1:kLim+1;
WMat(k+1,2:kLim+1) = cos(currExtrFreqs(1,k+1)'*([1:kLim]));
Num = ((-1).^k);
Denom = (valW4XFr(k+1));
WMat(k+1,kLim+2) = (Num./Denom)';
% end

AlphaDeltaVec = pinv(WMat)*(HdrVec)';
delta = AlphaDeltaVec(length(AlphaDeltaVec))

P = 0;
for pCtr = 1:length(AlphaDeltaVec)-1
    P = P + AlphaDeltaVec(pCtr)*(cos(FrGrid*(pCtr-1)) );
end

E = W.*([Hdr - P]);

if PassbandType==1 % lowpass
indicesLoBandOnFineGrid = 1:1:round(wp*LenGrid);
indicesHiBandOnFineGrid = round(ws*LenGrid):1:LenGrid;
elseif PassbandType==2
indicesLoBandOnFineGrid = 1:1:round(ws*LenGrid);
indicesHiBandOnFineGrid = round(wp*LenGrid):1:LenGrid;    
end

ftsz1 = 15;
ftsz2 = 13;

figure(56)
clf
subplot(311)
yarg = E(1,[indicesLoBandOnFineGrid]);
plot( (indicesLoBandOnFineGrid-1)/LenGrid, yarg )
ylabel('Wtd Error')
xlabel(['(a) Normalized Frequency, Band 1'])
line([0,(max(indicesLoBandOnFineGrid)-1)/LenGrid],[abs(delta), abs(delta)]);
line([0,(max(indicesLoBandOnFineGrid)-1)/LenGrid],[-abs(delta), -abs(delta)]);
plotmin = 1.2*min(yarg);
plotmax = 1.2*max(yarg);
axis([0,leftXitionBndLim,plotmin,plotmax])

subplot(312)
yarg = E(1,[indicesHiBandOnFineGrid]);
plot((indicesHiBandOnFineGrid-1)/LenGrid, yarg )
xlabel(['(b) Normalized Frequency, Band 2'])
ylabel('Wtd Error')
line([(min(indicesHiBandOnFineGrid)-1)/LenGrid,1],[abs(delta), abs(delta)]);
line([(min(indicesHiBandOnFineGrid)-1)/LenGrid,1],[-abs(delta), -abs(delta)]);
plotmin = 1.2*min(yarg);
plotmax = 1.2*max(yarg);
axis([rightXitionBndlim,1,plotmin,plotmax])

subplot(313)
fr = 20*log10( abs(Q.*P) + eps ) ;
plot([0:1:length(P)-1]/length(P),fr)
xlabel(['(c) Norm Freq, (L = ',num2str(L),')'])
ylabel('Mag, dB')
axis([0,1,max([min(fr),-130 ]),5])


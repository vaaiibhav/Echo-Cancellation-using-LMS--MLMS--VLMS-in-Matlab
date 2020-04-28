function DispPZCoordZeros

global hpzcoord1115928 hndUCAx115928  hndOneTwoOrFour hndPlotPZ1 hndPlotPZ2
global SR115928  hndPlotPZ1 hndPlotPZ2 hndPlotZero1 hndTxt1115928
global hndPlotPZ3 hndPlotPZ4
global dataPlotAxis1115928 dataPlotAxis2115928 dataPlotAxis3115928 dataPlotAxis4115928
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

ftszPole = 14;
ftszZero = 8;

ftsz1 = 14;
ftsz2 = 14;

currPt = get(hndUCAx115928,'CurrentPoint');

hold on
%scnsize = get(hndFig115928,'Position');

subplot(hndUCAx115928)

mag = sqrt(currPt(1,1)^2 + currPt(1,2)^2);

if mag > 1.1
   set(hpzcoord1115928,'String','Mag > 1.1 !');
   return
end

if abs(currPt(1,2)) < 0.005 % snap to real axis
   currPt(1,2) = 0;
end
   
ang = atan2(currPt(1,2),currPt(1,1))/pi;

stLoc = ['  Mag = ',num2str(mag,2),'; Norm Freq = ',num2str(ang,2)];

realPZ = currPt(1,1);
imagPZ = currPt(1,2);

if abs(realPZ)<0.005
    realPZ = 0;
end
if abs(imagPZ)<0.005
    imagPZ = 0;
end

if abs(1-realPZ)<0.005
    realPZ = 1;
end

if abs(1 - imagPZ)<0.005
    imagPZ = 1;
end

if abs(-1-realPZ)<0.005
    realPZ = -1;
end

if abs(-1-imagPZ)<0.005
    imagPZ = -1;
end

set(hpzcoord1115928,'String',stLoc);

IsOneTwoOrFour = get(hndOneTwoOrFour,'value');

SR115928 = 256;
Imp115928 = [1  zeros(1,255)];
    
if IsOneTwoOrFour==1  % one zero
      PZ(1) = realPZ + j*imagPZ; 
    
        set(hndPlotPZ1,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ2,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ3,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ4,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)     
        set(hndPlotZero1,'xdata',0,'ydata',0,'marker','x','markersize',ftszPole)
        set(hndTxt1115928,'string','1')
        a = 1;
        b = [1  -PZ(1)];  % a "zero" on the plot is the frequency at which the transfer fcn is to go to 
        % zero; to convert that frequency into a transfer function
        % coefficient, a negative sign is needed since 1 + b1*z^(-1)
        % implies that z^1 + b1  = 0 which implies that z = -b1 which
        % implies that b1 = -z (z being a frequency designated on the plot
        % at which the transfer function is to go to 0)       
        theImp115928Resp = filter(b,a,Imp115928);
        zarg = -pi:2*pi/255:pi;
            z = exp(j*zarg);
            theFreqResp = 1 -PZ(1)*z.^(-1);
        strzXform = ['1 - (',num2str(PZ(1),3),')z^{-1}'];
elseif IsOneTwoOrFour==2
    PZ(1) = realPZ + j*imagPZ;
    PZ(2) = realPZ - j*imagPZ;

        set(hndPlotPZ1,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ2,'xdata',realPZ,'ydata',-imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ3,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ4,'xdata',realPZ,'ydata',-imagPZ,'marker','o','markersize',ftszZero)    
        set(hndPlotZero1,'xdata',0,'ydata',0,'marker','x','markersize',ftszPole)
        set(hndTxt1115928,'string','2')
        a = 1;
        b = [1  -PZ(1)]; 
        theImp115928Resp = filter(b,a,Imp115928);
        b = [1  -PZ(2)];
        theImp115928Resp = filter(b,a,theImp115928Resp);
            bb = conv([1 -PZ(1)],[1  -PZ(2)]);
            zarg = -pi:2*pi/255:pi;
            z = exp(j*zarg);
            theFreqResp = bb(1) + bb(2)*z.^(-1)  + bb(3)*z.^(-2);
        lowvals = find(abs(bb)<0.005);
        bb(lowvals) = 0;
        
strzXform = ['1 + (',num2str(bb(2),3),')z^{-1}',' + (',num2str(bb(3),3),')z^{-2}'];
else  % number of zeros = 4
    PZ(1) = realPZ + j*imagPZ;
    PZ(2) = realPZ - j*imagPZ;
    newmag = 1/sqrt(realPZ^2 + imagPZ^2);
    PZ(3) = PZ(1)*(newmag^2);
    PZ(4) = PZ(2)*(newmag^2);
    
  set(hndPlotPZ1,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ2,'xdata',realPZ,'ydata',-imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotZero1,'xdata',0,'ydata',0,'marker','x','markersize',ftszPole)
        set(hndPlotPZ3,'xdata',real(PZ(3)),'ydata',imag(PZ(3)),'marker','o','markersize',ftszZero)
        set(hndPlotPZ4,'xdata',real(PZ(4)),'ydata',imag(PZ(4)),'marker','o','markersize',ftszZero)      
        set(hndTxt1115928,'string','4')
a = 1;
bb = poly([PZ(1),PZ(2),PZ(3),PZ(4)]);          
theImp115928Resp = filter(bb,a,Imp115928);
zarg = -pi:2*pi/255:pi;
z = exp(j*zarg);
theFreqResp = bb(1) + bb(2)*z.^(-1) + bb(3)*z.^(-2) + bb(4)*z.^(-3) + bb(5)*z.^(-4);
lowvals = find(abs(bb)<0.005);
bb(lowvals) = 0;
bb = real(bb);     
strzXform = ['1 + (',num2str(bb(2),3),')z^{-1}',' + (',num2str(bb(3),3),')z^{-2}',' + (',num2str(bb(4),3),')z^{-3}',' + z^{-4}'];
end

plot1max = 1.2*max(abs(theFreqResp));

subplot(dataPlotAxis1115928)
cla
hold on
x = (1/128)*(-127:1:128);
plot(x,abs(theFreqResp),'k')
set(gca,'fontsize',ftsz1)
xlabel('(a) Normalized Freq','fontsize',ftsz2)
ylabel('Mag','fontsize',ftsz2)
axis([-1  1  0  plot1max])

theangle = unwrap(angle(theFreqResp));
plot2max = 1.2*max(abs(theangle));

subplot(dataPlotAxis2115928)
cla

plot(x,theangle,'k')
set(gca,'fontsize',ftsz1)
xlabel('(b) Normalized Freq','fontsize',ftsz2)
ylabel('Rad','fontsize',ftsz2)
axis([-1  1  min([-0.5  -plot2max])  max([0.5   plot2max])    ])

plot3max = 1.2*max(abs(real(theImp115928Resp(1,1:32))));

subplot(dataPlotAxis3115928)
cla

set(gca,'fontsize',ftsz1)
stem(real(theImp115928Resp),'k.')
ylabel('Real','fontsize',ftsz2)
xlabel('(c) Sample, Impulse Response','fontsize',ftsz2)
axis([0  32   min([-1  -plot3max])  max([1  plot3max])   ])

plot4max = 1.2*max(abs(imag(theImp115928Resp(1,1:32))));
subplot(dataPlotAxis4115928)
cla

set(gca,'fontsize',ftsz1)
stem(imag(theImp115928Resp),'k.')
ylabel('Imag','fontsize',ftsz2)
xlabel('(d) Sample, Impulse Response','fontsize',ftsz2)
axis([0  32   min([-1  -plot4max])  max([1  plot4max])  ])

subplot(hndUCAx115928)
hndT = title(strzXform);
set(hndT,'fontsize',12)

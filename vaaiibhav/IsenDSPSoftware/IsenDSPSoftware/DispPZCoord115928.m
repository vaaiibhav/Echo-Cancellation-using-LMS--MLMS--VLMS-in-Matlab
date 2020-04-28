function DispPZCoord115928

global hpzcoord1115928 hndUCAx115928  hndPoleOrZero115928 hndPlotPZ1 hndPlotPZ2
global dataPlotAxis1115928 dataPlotAxis2115928 dataPlotAxis3115928 dataPlotAxis4115928
global SR115928 hndOneOrTwo115928 hndPlotPZ1 hndPlotPZ2 hndPlotZero1 hndPlotZero2 hndTxt1115928
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

ftsz2 = 15;
ftsz1 = 13;
ftszPole = 14;
ftszZero = 8;

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

% do plots first
if(get(hndOneOrTwo115928,'value')==0)
    lenPZ=1;
    PZ(1) = realPZ + j*imagPZ; 
else
    lenPZ=2;
    PZ(1) = realPZ + j*imagPZ;
    PZ(2) = realPZ - j*imagPZ;
end
 
isPorZ = get(hndPoleOrZero115928,'value');

SR115928 = 256;
Imp115928 = [1  zeros(1,255)];
    
if isPorZ==1  % pole
    if lenPZ==1
       set(hndPlotPZ1,'xdata',realPZ,'ydata',imagPZ,'marker','x','markersize',ftszPole)
       set(hndPlotPZ2,'xdata',realPZ,'ydata',imagPZ,'marker','x','markersize',ftszPole)
       
       set(hndPlotZero1,'xdata',0,'ydata',0,'marker','o','markersize',ftszZero)
       set(hndPlotZero2,'xdata',0,'ydata',0,'marker','o','markersize',ftszZero)
        set(hndTxt1115928,'string','1')
        b = 1;
        a = [1  -PZ(1)];
        theImp115928Resp = filter(b,a,Imp115928);
        if (abs(PZ(1))>=1.0)
            theFreqResp = fftshift(fft(theImp115928Resp));
        else
            zarg = -pi:2*pi/255:pi;
            z = exp(j*zarg);
            theFreqResp = 1./(1 -PZ(1)*z.^(-1));  
        end
        strzXform = ['1/(1 - (',num2str(PZ(1),3),')z^{-1}'];
    else
       set(hndPlotPZ1,'xdata',realPZ,'ydata',imagPZ,'marker','x','markersize',ftszPole)
       set(hndPlotPZ2,'xdata',realPZ,'ydata',-imagPZ,'marker','x','markersize',ftszPole) 
       
       set(hndPlotZero1,'xdata',0,'ydata',0,'marker','o','markersize',ftszZero)
       set(hndPlotZero2,'xdata',0,'ydata',0,'marker','o','markersize',ftszZero)
        set(hndTxt1115928,'string','2')
        b = 1;
        a = [1  -PZ(1)];      
        a2 = [1  -PZ(2)];
        theImp115928Resp = filter(b,a,Imp115928);
        theImp115928Resp = filter(b,a2,theImp115928Resp);
        if (abs(PZ(1))>=1.0||abs(PZ(2))>=1.0) 
        theFreqResp = fftshift(fft(theImp115928Resp));
        else % use z-transform to get freq response
            aa = conv([1 -PZ(1)],[1  -PZ(2)]);
            zarg = -pi:2*pi/255:pi;
            z = exp(j*zarg);
            theFreqResp = 1./(aa(1) + aa(2)*z.^(-1)  + aa(3)*z.^(-2));   
        end
        aa = conv([1 -PZ(1)],[1  -PZ(2)]);
        lowvals = find(abs(aa)<0.005);
        aa(lowvals) = 0;
        strzXform = ['1/(1 + (',num2str(aa(2),3),')z^{-1}',' + (',num2str(aa(3),2),')z^{-2}'];
     end
else  % zero
    if lenPZ==1

        set(hndPlotPZ1,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ2,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
       
        set(hndPlotZero1,'xdata',0,'ydata',0,'marker','x','markersize',ftszPole)
        set(hndPlotZero2,'xdata',0,'ydata',0,'marker','x','markersize',ftszPole)
        set(hndTxt1115928,'string','1')
        a = 1;
        b = [1  -PZ(1)];  % a "zero" on the plot is the frequency at which the transfer fcn is to go to 
        % zero; to convert that frequency into a transfer function
        % coefficient, a negative sign is needed since 1 + b1*z^(-1)
        % implies that z^1 + b1  = 0 which implies that z = -b1 which
        % implies that b1 = -z (z being a frequency designated on the plot
        % at which the transfer function is to go to 0)       
        theImp115928Resp = filter(b,a,Imp115928);
        %theFreqResp = fftshift(fft(theImp115928Resp));
        zarg = -pi:2*pi/255:pi;
            z = exp(j*zarg);
            theFreqResp = 1 -PZ(1)*z.^(-1);
        strzXform = ['1 - (',num2str(PZ(1),3),')z^{-1}'];
    else

        set(hndPlotPZ1,'xdata',realPZ,'ydata',imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotPZ2,'xdata',realPZ,'ydata',-imagPZ,'marker','o','markersize',ftszZero)
        set(hndPlotZero1,'xdata',0,'ydata',0,'marker','x','markersize',ftszPole)
        set(hndPlotZero2,'xdata',0,'ydata',0,'marker','x','markersize',ftszPole)
        set(hndTxt1115928,'string','2')
        a = 1;
        b = [1  -PZ(1)]; 
        theImp115928Resp = filter(b,a,Imp115928);
        b = [1  -PZ(2)];
        theImp115928Resp = filter(b,a,theImp115928Resp);
        % theFreqResp = fftshift(fft(theImp115928Resp));
            bb = conv([1 -PZ(1)],[1  -PZ(2)]);
            zarg = -pi:2*pi/255:pi;
            z = exp(j*zarg);
            theFreqResp = bb(1) + bb(2)*z.^(-1)  + bb(3)*z.^(-2);
        lowvals = find(abs(bb)<0.005);
        bb(lowvals) = 0;
        
        strzXform = ['1 + (',num2str(bb(2),3),')z^{-1}',' + (',num2str(bb(3),2),')z^{-2}'];
    end
end
plot1max = 1.2*max(abs(theFreqResp));

subplot(dataPlotAxis1115928)
cla
hold on
x = (1/128)*(-127:1:128);
plot(x,abs(theFreqResp))
set(gca,'fontsize',ftsz1)
xlabel('(a) Normalized Freq','fontsize',ftsz2)
ylabel('Mag','fontsize',ftsz2)
axis([-1  1  0  plot1max])

theangle = unwrap(angle(theFreqResp));
plot2max = 1.2*max(abs(theangle));

subplot(dataPlotAxis2115928)
cla

plot(x,theangle)
set(gca,'fontsize',ftsz1)
xlabel('(b) Normalized Freq','fontsize',ftsz2)
ylabel('Rad','fontsize',ftsz2)
axis([-1  1  min([-0.5  -plot2max])  max([0.5   plot2max])    ])

plot3max = 1.2*max(abs(real(theImp115928Resp(1,1:32))));

subplot(dataPlotAxis3115928)
cla

set(gca,'fontsize',ftsz1)
stem(real(theImp115928Resp),'b.')
ylabel('Real','fontsize',ftsz2)
xlabel('(c) Impulse Response','fontsize',ftsz2)
axis([0  32   min([-1  -plot3max])  max([1  plot3max])   ])

plot4max = 1.2*max(abs(imag(theImp115928Resp(1,1:32))));
subplot(dataPlotAxis4115928)
cla

set(gca,'fontsize',ftsz1)
stem(imag(theImp115928Resp),'b.')
ylabel('Imag','fontsize',ftsz2)
xlabel('(d) Impulse Response','fontsize',ftsz2)
axis([0  32   min([-1  -plot4max])  max([1  plot4max])  ])

subplot(hndUCAx115928)
hndT = title(texlabel(strzXform));
set(hndT,'fontsize',13)

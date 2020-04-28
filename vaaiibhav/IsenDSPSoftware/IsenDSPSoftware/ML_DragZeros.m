function ML_DragZeros
%function ML_DragZeros
%
%Opens a GUI having a Unit Circle and various axes. Mouse can be used to
%move cursor in vicinity of Unit Circle and one or two or four zeros
%(user-selectable number) will be created and will move with the cursor. The
%corresponding impulse response, frequency response, and z-transform will be
%displayed. 
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

global hpzcoord1115928    
global dataPlotAxis1115928 dataPlotAxis2115928 dataPlotAxis3115928 dataPlotAxis4115928 
global hndUC115928 hndCenUC115928 hndRealAxUC115928 hndImagAxUC115928 
global hndUCLabel115928 hndUCAx115928  Imp115928 
global hndFig115928 hndOneTwoOrFour SR115928 hndPlotPZ3 hndPlotPZ4
global hndPlotPZ1 hndPlotPZ2 hndPlotZero1 hndTxt1115928

SR115928 = 256;
Imp115928 = [1  zeros(1,31)];

scnsize = get(0,'ScreenSize');
pos1 = [0.01*scnsize(3), 0.05*scnsize(4), 0.98*scnsize(3), 0.86*scnsize(4)];
 
hndFig115928 = figure(115928);
set(hndFig115928,'Color',[1 1 1],'DoubleBuffer','On','paperpositionmode','auto','position',pos1)
set(hndFig115928,'numbertitle','off','name','Z-Plane Plot with Corresponding Impulse & Chirp Responses')
clf

SR = 256;
theta = 0:(2*pi)/SR:(2*pi-(2*pi)/SR);

xs = cos(theta);
ys = sin(theta);

horlinex = [1,(-1.2:0.1:1.2)];

vertliney = zeros(size(horlinex));

Vliney = -1.1:0.1:1.1;
Vlinex = zeros(size(Vliney));

hndUCAx115928 = axes('position',[0.525 0.15 0.45 0.6]);

hold on
hndUC115928 = plot(xs,ys,'k:');

hndCenUC115928 = plot(0,0,'k.');

hndRealAxUC115928 = plot(horlinex,vertliney,'k:');
hndImagAxUC115928 = plot(Vlinex, Vliney,'k:');
set(gca,'fontsize',13)
ylabel('Imaginary','fontsize',13)
xlabel('(e) Real','fontsize',13)

hndUCLabel115928 = text(0.5,-0.97,'Unit Circle','fontsize',13);
axis([-1.6, 1.6, -1.5, 1.5])

scnsz = get(0,'ScreenSize');
set(hndFig115928,'Position',[0.005*scnsz(3) 0.04*scnsz(4) 0.99*scnsz(3) 0.9*scnsz(4)]);
width = 0.95*scnsz(3);

posSetPoleZero115928 = [0.5*width 0.0475*scnsz(4) 0.1*width 0.035*scnsz(4)];
poslabelCheckbox = [0.68*width 0.0125*scnsz(4) 0.15*width 0.035*scnsz(4)];
% uicontrol('style','text','position',poslabelCheckbox,'string','Select No. of Zeros','BackgroundColor',[1 1 1]);
hndOneTwoOrFour = uicontrol(hndFig115928,'Style','popup','Position',posSetPoleZero115928,'string','Single|CC Pair|Quad','backgroundcolor',[1 1 1]); 

hpzcoord1115928 = uicontrol(hndFig115928,'style','text','position',...
  [0.54*scnsz(3) 0.795*scnsz(4) 0.4*scnsz(3) 0.04*scnsz(4)]);
set(hpzcoord1115928,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','fontsize',13);

dataPlotAxis1115928 = axes('Position',[0.12 0.8 0.32 0.12]);

dataPlotAxis2115928 = axes('Position',[0.12 0.585 0.32 0.12]);

dataPlotAxis3115928 = axes('Position',[0.12 0.3675 0.32 0.12]);

dataPlotAxis4115928 = axes('Position',[0.12 0.15 0.32 0.12]);

figure(hndFig115928)
subplot(hndUCAx115928)

realPZ = 0;
imagPZ = 0;

hndPlotPZ1 = plot(realPZ,imagPZ,'ko','markersize',8);
hndPlotPZ2 = plot(realPZ,imagPZ,'ko','markersize',8);
hndPlotPZ3 = plot(realPZ,imagPZ,'ko','markersize',8);
hndPlotPZ4 = plot(realPZ,imagPZ,'ko','markersize',8);

hndPlotZero1 = plot(realPZ,imagPZ,'kx','markersize',8);
% hndPlotZero2 = plot(realPZ,imagPZ,'kx','markersize',8);
hndTxt1115928 = text(0.06,0.06,'');
set(hndFig115928,'WindowButtonMotionFcn','DispPZCoordZeros')





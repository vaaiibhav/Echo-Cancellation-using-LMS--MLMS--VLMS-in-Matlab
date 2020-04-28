function ML_Correlation
%function ML_Correlation
%
%This program demonstrates how correlation sequences are
%computed, one sample at a time. One function "slides" over the
%other, one sample at a time. When one function is centered over      
%the other, the offset is zero; the number of samples of offset
%is referred to as the lag number or index.               
%At each lag, the product of all overlapping samples is taken 
%and all products are summed to yield the correlation value
%for the particular lag.  
%The plot of such numbers is known as the correlation sequence.                  
%Each time a computation is performed, the program pauses   
%for you to observe what has happened.  To view the next
%computational result in the sequence, press any key. To exit any 
%MATLAB program, including this one, at any time, press
%Ctrl-C.  Remember to press any key to perform the next computation, and press any 
%key twice after each of the first two correlation examples are finished to move to
%the next example (there are three separate examples executed in sequence).
%
% The call:
%
% ML_Correlation
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

%-------------AutoCorrelation of a Rectangle---------------------------
yshift(1,1:1:23) = zeros;
CorrelSeqVal(1,1:15) = zeros;
x(1,1:23) = zeros;
y(1,1:40) = zeros;
x(1,8:15) = 0.8*[1 1 1 1 1 1 1 1];
y(1,16:23)= 1.25*[1 1 1 1 1 1 1 1];

for a = 0:1:16  
   
yshift(1,1:23) = y(1,a+1:23+a);
CorrelSeqVal(1,a+1) = sum(x.*yshift);

scnsize = get(0,'ScreenSize');
pos2 = [0.0025*scnsize(3), 0.06*scnsize(4),...
      0.995*scnsize(3), 0.84*scnsize(4)];
figure(196);
clf
set(196,'Name','Correlation Sequence of Two Rectangles','Color',[1 1 1],...
   'Position',pos2,'NumberTitle','off', 'DoubleBuffer','on');
subplot(2,1,1)

hold off
hndPlot = stem(x,'ro');
set(hndPlot,'markersize',5)
hold on
hndPlot = stem(yshift,'bo');
set(hndPlot,'markersize',5)

set(gca,'fontsize',18)
ylabel(['Amplitude'],'fontsize',18)
xlabel(['(a)  Sample Index'],'fontsize',18)
axis([0 25 0 1.8])

subplot(2,1,2)
hndPlot = stem(-8:-8+a,CorrelSeqVal(1,1:a+1),'bo');
set(hndPlot,'markersize',5)
axis([-8 8 -15 15])
if abs(CorrelSeqVal(1,a+1))<10^(-6)
    CorrelSeqVal(1,a+1) = 0;
end
text(-5.5,11,['Correlation Sequence Value at Lag No. ',num2str(a-8),' is '...
      ,num2str(CorrelSeqVal(1,a+1),3)],'fontsize',16)

set(gca,'fontsize',18)
ylabel(['Amplitude'],'fontsize',18)
xlabel(['(b)  Lag Number (Press Any Key for Next Lag)'],'fontsize',18)
if a == 16
xlabel(['(b)  Lag Number (Press Any Key for Next Demo)'],'fontsize',18)    
end

pause
end
xlabel(['(b)  Lag Number (Press Any Key for Next Demo)'],'fontsize',18)
pause
%-------------------AutoCorrelation of a Sinewave----------------------------------
set(196,'Name','Autocorrelation Seq. of a Sinewave')
yshift(1,1:1:23) = zeros;
CorrelSeqVal(1,1:15) = zeros;
x(1,1:23) = zeros;
y(1,1:40) = zeros;
t=0:1:7;
x(1,8:15) = sin(2*pi*t/8);
y(1,16:23)= sin(2*pi*t/8);

for a = 0:1:16
yshift(1,1:23) = y(1,a+1:23+a);
CorrelSeqVal(1,a+1) = sum(x.*yshift);

subplot(2,1,1)
hold off
plot(x,'r')
hold on
hndPlot = stem(x,'ro');
set(hndPlot,'markersize',5)
hndPlot = stem(yshift,'bd');
set(hndPlot,'markersize',6)

set(gca,'fontsize',18)
ylabel(['Amplitude'],'fontsize',18)
xlabel(['(a)  Sample Index'],'fontsize',18)
axis([0 25 -1.8 1.8])

subplot(2,1,2)
hndPlot = stem(-8:-8+a,CorrelSeqVal(1,1:a+1),'bo');
set(hndPlot,'markersize',5)
axis([-8 8 -15 15])
if abs(CorrelSeqVal(1,a+1))<10^(-6)
    CorrelSeqVal(1,a+1) = 0;
end
text(-5.5,11,['Correlation Sequence Value at Lag No. ',num2str(a-8),' is '...
      ,num2str(CorrelSeqVal(1,a+1),3)],'fontsize',16)

set(gca,'fontsize',18)
ylabel(['Amplitude'],'fontsize',18)
xlabel(['(b)  Lag Number (Press Any Key for Next Lag)'],'fontsize',18)
if a == 16
xlabel(['(b)  Lag Number (Press Any Key for Next Demo)'],'fontsize',18)    
end
pause
end
xlabel(['(b)  Lag Number (Press Any Key for Next Demo)'],'fontsize',18)
pause

%---------------Correlation of a Sine and Cosine-------------------
set(196,'Name','Correlation Seq. of Sine & Cosine')
yshift(1,1:1:23) = zeros;
CorrelSeqVal(1,1:15) = zeros;
x(1,1:23) = zeros;
y(1,1:40) = zeros;
t=0:1:7;
x(1,8:15) = sin(2*pi*t/8);
y(1,16:23)= cos(2*pi*t/8);

for a = 0:1:16
yshift(1,1:23) = y(1,a+1:23+a);
CorrelSeqVal(1,a+1) = sum(x.*yshift);

subplot(2,1,1)
hold off
plot(x,'r')
hold on
hndPlot = stem(x,'ro');
set(hndPlot,'markersize',5)
hndPlot = stem(yshift,'bd');
set(hndPlot,'markersize',6)
axis([0 25 -1.8 1.8])
set(gca,'fontsize',18)
ylabel(['Amplitude'],'fontsize',18)
xlabel(['(a)  Sample Index'],'fontsize',18)

subplot(2,1,2)
hndPlot = stem(-8:-8+a,CorrelSeqVal(1,1:a+1),'bo');
set(hndPlot,'markersize',4)
axis([-8 8 -15 15])
if abs(CorrelSeqVal(1,a+1))<10^(-6)
    CorrelSeqVal(1,a+1) = 0;
end
text(-5.5,11,['Correlation Sequence Value at Lag No. ',num2str(a-8),...
   ' is  ',num2str(CorrelSeqVal(1,a+1),3)],'fontsize',16)
set(gca,'fontsize',18)
ylabel(['Amplitude'],'fontsize',18)
xlabel(['(b)  Lag Number (Press Any Key for Next Lag)'],'fontsize',18)
if a == 16
xlabel(['(b)  Lag Number (Demo Finished!)'],'fontsize',18)
return
end
pause
end


clc;
clear all;
%warning off;
close all;
 


load farspeech
ss=x;
lengthx = abs(real(1:length(x)));
%ss= audioread('C:\Users\OM\Desktop\clean\sp12.wav');
xlabel('Time [sec]');
ylabel('Amplitude');
title('Far-End Speech Signal');
figure(1)
%subplot (3,3,1);
plot(x);
wavwrite(x,'farspeech');

load nearspeech;
xlabel('Time [sec]');
ylabel('Amplitude');
title('Near-End Speech Signal');
figure(3)
%subplot (3,3,3);
plot(v);
wavwrite(v,'nearend');
lengthv = abs(real(1:length(v)));
x=x(1:lengthv);
v=v(1:lengthx);

NrEndPlusFarEnd(n) = v(1,lengthv) + x(1,lengthx);
n=audioread('echoedsignal.wav');
xlabel('Time [sec]');
ylabel('Amplitude');
title('Near end and Far-End Speech Signal');
figure(2)
plot(NrEndPlusFarEnd);% subplot (1,2,2), disp (n);


xlabel('Time [sec]');
ylabel('Amplitude');
title('Microphone Speech Signal');
%subplot (3,3,4);
figure(4)
plot(a);
wavwrite(a,'microphone');

% plot(a);
%Initialization
 N=20115;
 %Hpsd=dspdata.psd(N);
p=1024;
 mu=0.04;
w=zeros(p,1);
% x=zeros(N,1);
% d=zeros(N,1);
% Input signals 
W0 = zeros(1,2048);
x = x(1:length(W0)*floor(length(x)/length(W0)));
d = a(1:length(W0)*floor(length(a)/length(W0)));

 
%  figure;
%  plot(s);
%Algorithm
for i=p:N
    xvec=n(i:-1:i-p+1);
    y(i)=w'*xvec;
     e(i)=a(i)-y(i);
    w=w+mu*e(i)*xvec;
end

[c2,e2] = filter(e,x,d)
wavwrite(e2,'error_signal');
%plot e 
xlabel('Time [sec]');
ylabel('Amplitude');
title('error_signal Signal');
%subplot (3,3,5);
figure(6)
plot(e2);
% title('');
%wavwrite(MSE,'MSE_Signal');

% plot((s'-e));
% xlabel('time index'); 
% ylabel('signal value');
%Calculating MSE
[R C]=size(ss);
newmatrix(i)= ss(i)-a(i);
newmatrix2(i) = ss(i)-e(i);
for i=1:N
    err(i) = (((ss(i)-e(i)).^2)/(R*C));
%      err(i)=(ss(i)-e(i)).^2;
    % nn(i)=n(i).^2;
%     err(i)=(s(i)-e(i));
%      nn(i)=n(i);
%     ss(i)=s(i);
end
mse =psnr1(ss);

 % size of signals must be same length 
MSE=sqrt(err);
% disp (MSE);

%wavwrite(err,'error_signal');

% MSE=mean(err)
xlabel('Time [sec]');
ylabel('Amplitude');
title('err Signal');
%subplot (3,3,5);
figure(5)
plot(MSE);
title('mse');
%wavwrite(MSE,'MSE_Signal');

% %  Calculating SNR INPUT
% SNR=snr(e,(ss-e))
 rms_signal=sqrt(mean(a.^2)).'
%    z = ss .* a.';
%   reshape(ss,size(a));
%   reshape(a,size(ss));
  rms_echo=sqrt(mean(newmatrix.^2));
  Lsig=10*log10(rms_signal);
  Lech=10*log10(rms_echo);
  ERLEo=Lsig-Lech
%   xlabel('n');
% ylabel('ERLE');
% title('ERLE Signal');
% %subplot (3,3,4);
% figure(6)
% plot(ERLEo);
%   wavwrite(ERLEo, 'ERLE_SIGNAL');

% 
% %  
% % Calculating SNR OUTPUT
% % SNR=snr(e,(ss-e))
 rms_signal=sqrt(mean(e.^2));
  rms_echo=sqrt(mean(newmatrix2.^2));
  Lsig=10*log10(rms_signal);
  Lech=10*log10(rms_echo);
  ERLEi=Lsig-Lech

function [Type]= LVxAmp_V_MagResponse(Imp,FreqRange,IncMag,LogPlot)
% function [Type]= LVxAmp_V_MagResponse(Imp,FreqRange,IncMag,LogPlot)
% Imp is a linear phase impulse response of Types I, II, III, or IV
% Pass FreqRange as 0 for -pi to pi; 1 for -2pi to 2pi; 2 for 0 to pi; 3
% for 0 to 2pi, and 4 for 0 to 4pi
% Pass IncMag as 1 to include magnitude plots, or 0 for amplitude plots
% only
% Pass LogPlot as 0 for linear magnitude plot or 1 for 20log10(Mag) plot;
% The output variable Type is returned as 1,2,3, or 4 for Types I, II, III,
% or IV, respectively, or 0 if Imp is not a linear phase impulse response.
% Test calls:
% [Type]= LVxAmp_V_MagResponse([1,0,1],0,1,0)
% [Type]= LVxAmp_V_MagResponse([1,1,1,1],0,1,0)
% [Type]= LVxAmp_V_MagResponse([1,0,-1],0,1,0)
% [Type]= LVxAmp_V_MagResponse([1,1,-1,-1],0,1,0)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

% identify which Type
Type = 0;
L = length(Imp);
M = (L-1)/2;

T12test = (sum((Imp - fliplr(Imp)).^2 ))^0.5;
T34test = (sum((Imp + fliplr(Imp)).^2 ))^0.5;

if T34test < 10^(-15) % Type III or IV
    if rem(L,2)==0 % even length
        Type = 4
    else % odd length
        Type = 3
    end
elseif T12test < 10^(-15) % it's Type I or II
        if rem(L,2)==0 % even length
            Type = 2
        else
            Type = 1
        end
else % not any of the four Types
    Message = 'Impulse response is not Type I, II, III, or IV. Exiting procedure.'
    Type = 0
    return
end

ftsz1 = 16;
ftsz2 = 14;

scnsize = get(0,'ScreenSize');
pos1 = [0.005*scnsize(3), 0.04*scnsize(4), 0.99*scnsize(3), 0.86*scnsize(4)];

figure(655)
set(655,'Color',[1 1 1],'Position',pos1,'Name',['Amplitude V. Magnitude Respones for Linear Phase Filters'],'Numbertitle','off')

clf

incF = 2*pi/512; 

if FreqRange==0
    argF = -pi+incF:incF:pi;
elseif FreqRange==1
    argF = -2*pi:incF:2*pi;
elseif FreqRange==2
    argF = 0:incF:pi;
elseif FreqRange==3
    argF = 0:incF:2*pi;
elseif FreqRange==4
    argF = 0:incF:4*pi;
end

if Type == 1
    Hr = 0;
    for n = 1:1:M
        Hr = Hr + 2*Imp(M-n+1)*cos(argF*n);
    end
    Hr = Hr + Imp(M+1);
    Hph = -(argF)*M;

elseif Type==2
    Hr = 0;
    %for n = 1:1:L/2
    %    Hr = Hr + 2*Imp(L/2 - n + 1)*cos(argF*(n - 0.5));
    %end
    for n = 0:1:L/2-1
        Hr = Hr + 2*Imp(n+1)*cos(argF*(M-n));
    end
    Hph = -argF*M;

elseif Type==3
    Hr = 0;
    for n = 0:1:M-1
        Hr = Hr + 2*Imp(n+1)*sin(argF*(M-n));
    end
Hph = pi/2 - argF*M;
 
elseif Type==4
        Hr = 0;
    for n = 0:1:L/2-1
        Hr = Hr + 2*Imp(n+1)*sin(argF*(M-n));
    end
Hph = pi/2 - argF*M;

end

if IncMag==1
    subplot(221); 
    plot(argF/pi, Hr)
    set(gca,'fontsize',ftsz2)
    xlabel('(a) Normalized Frequency','fontsize',ftsz1)
    ylabel('Amplitude','fontsize',ftsz1)
    
    subplot(222); 
    plot(argF/pi,Hph)
    set(gca,'fontsize',ftsz2)
    xlabel('(b) Normalized Frequency','fontsize',ftsz1)
    ylabel('Radians','fontsize',ftsz1)
    
    z = exp(-j*argF);
    Hz = polyval(Imp,z)./polyval(1,z); 
     
    subplot(223); 
    
    if LogPlot==1
    plot(argF/pi, 20*log10(abs(Hz)+10^(-5)) );
        ylabel('Mag (dB)','fontsize',ftsz1)
    else
    plot(argF/pi, abs(Hz));
        ylabel('Mag','fontsize',ftsz1)
    end

    set(gca,'fontsize',ftsz2)
    xlabel('(c) Normalized Frequency','fontsize',ftsz1)

    
    subplot(224); 
    plot(argF/pi,angle(Hz))
    set(gca,'fontsize',ftsz2)
    xlabel('(d) Normalized Frequency','fontsize',ftsz1)
    ylabel('Radians','fontsize',ftsz1)
else
    subplot(211); 
    plot(argF/pi, Hr)
    set(gca,'fontsize',ftsz2)
    xlabel('(a) Normalized Frequency','fontsize',ftsz1)
    ylabel('Amplitude','fontsize',ftsz1)
    
    subplot(212); 
    plot(argF/pi,Hph)
    set(gca,'fontsize',ftsz2)
    xlabel('(b) Normalized Frequency','fontsize',ftsz1)
    ylabel('Radians','fontsize',ftsz1) 
end


   



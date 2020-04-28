function OutputWave = LVxEncodeDTMF(DialedDigits,NoiseVar)
% function OutputWave = LVxEncodeDTMF(DialedDigits,NoiseVar)
%
% DialedDigits must consist of either 
% 1) the numbers 0,1,2...9 (no alphabetic
% characters or other non-numeric characters), or
% 2) all string format, i.e., all numbers and letters enclosed in single
% parentheses, in the following sample call:
% OW = LVxEncodeDTMF(['5';'4';'0';'A';'C';'6';'3'],0.05);
% OW = LVxEncodeDTMF(['5','4','0','A','C','6','3'],0.05);
%Sample Calls:
%
% OW = LVxEncodeDTMF([7 0 3 ],0.05)  % numeric data only
% OW = LVxEncodeDTMF(['5';'4';'0';'A';'C';'6';'3'],0.05); % alphanumeric data
% OW = LVxEncodeDTMF(['5','4','0','A','C','6','3'],0.05); % alphanumeric data
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
RowMat = [697; 770;852;941];
ColMat = [1209; 1336;1477;1633];

OutputWave = [];
for ctr = 1:1:length(DialedDigits) % identify row and col for each digit to generate two freqs
    if (isnumeric(DialedDigits(ctr))==1)
        dd = num2str(DialedDigits(ctr));
    else 
        dd = DialedDigits(ctr);
    end

if strcmp(dd,'1')==1
    rc = [1,1];
elseif strcmp(dd,'2')==1
    rc = [1,2];
elseif strcmp(dd,'3')==1
    rc = [1,3];
elseif strcmp(dd,'4')==1
    rc = [2,1];
elseif strcmp(dd,'5')==1
    rc = [2,2];
elseif strcmp(dd,'6')==1
    rc = [2,3];
elseif strcmp(dd,'7')==1
    rc = [3,1];
elseif strcmp(dd,'8')==1
    rc = [3,2];
elseif strcmp(dd,'9')==1
    rc = [3,3];
elseif strcmp(dd,'0')==1
    rc = [4,2];
elseif strcmp(dd,'A')==1
    rc = [1,4];
elseif strcmp(dd,'B')==1
    rc = [2,4];
elseif strcmp(dd,'C')==1
    rc = [3,4];
elseif strcmp(dd,'D')==1
    rc = [4,4];
elseif strcmp(dd,'*')==1
    rc = [4,1];
elseif strcmp(dd,'#')==1
    rc = [4,4];
end
f1 = RowMat(rc(1));
f2 = ColMat(rc(2));
NewTone = cos(2*pi*(0:1:7999)*f1/8000) + cos(2*pi*(0:1:7999)*f2/8000);
OutputWave = [OutputWave NewTone(1,1:2000) zeros(1,1000)];
end
OutputWave = OutputWave + NoiseVar*randn(1,length(OutputWave));
% plot(OutputWave)
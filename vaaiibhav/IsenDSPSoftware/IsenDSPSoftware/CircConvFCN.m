function theCircCon = CircConvFCN(Seq1,Seq2)
% function theCircCon = CircConvFCN(Seq1,Seq2)
% Performs circular convolution on Seq1 and Seq2
% 
% Example Calls:
% CircConvFCN([1 1 1 0 0 0 0 0],[110 121 115 132 115 107 143 117])
% CircConvFCN([1 1 1 1 -1 -1 -1 -1 zeros(1,8)],[0.25 0.5 0.75 1 1 0.75 0.5 0.25 zeros(1,8)])
% CircConvFCN([1 1 1 1 -1 -1 -1 -1],[0.25 0.5 0.75 1 1 0.75 0.5 0.25])

% CircConvFCN([1 1 1 0 0 0 0 0],[110+i*105 121 115 132 115 107 i*143 117])
% CircConvFCN([1 1*i 1 0 0 0 0 0],[110+i*105 121 115 132 115 107 i*143 117])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~length(Seq1)==length(Seq2)
   a = 'Sequences Not Same Length!'
   return
end

SR = length(Seq1);
Seq2 = [Seq2(1,1:1),fliplr(Seq2(1,2:length(Seq2)))];

MCircConv = zeros(1,length(Seq1));
maxplotSeq1 = 1.2*max(abs(Seq1));
maxplotSeq2 = 1.2*max(abs(Seq2));
maxplotCircConv = 1.2*max(abs(conv(Seq1,Seq2)));
lagctr = (-length(Seq1)/2) + 1;

for ctr = 1:1:length(Seq2)  
MCircConv(1,ctr:ctr) = sum(Seq1.*(Seq2));

xxx = Seq2(1,length(Seq2):length(Seq2)); % rotates to the right, or counterclockwise
Seq2(1,2:length(Seq2)) = Seq2(1,1:length(Seq2)-1);
Seq2(1,1:1) = xxx;

lagctr =lagctr + 1;

end
theCircCon =  MCircConv;


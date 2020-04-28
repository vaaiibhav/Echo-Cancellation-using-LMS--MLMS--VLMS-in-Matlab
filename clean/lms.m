function hFDAF = lms( lmsa,lmsb,lmsc,lmsd,lmse )
%LMS_ECHO Summary of this function goes here
%   Detailed explanation goes here




 hFDAF = adaptfilt.fdaf(lmsa, lmsb, lmsc, lmsd, lmse);

end


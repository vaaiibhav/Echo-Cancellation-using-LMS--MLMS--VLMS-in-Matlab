function SuccessiveApproximationResize
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
global selStepTypeSAR SARCompute hndWFSel hndBiasSet 
global hndLSBAmpSet hndSetPeakValTestSig hndFigSuccApprox hndSRSet

figpos = get(hndFigSuccApprox,'Position');
figwidth = figpos(3);
figheight = figpos(4);

 posSetPeakValTestSig = [0.08*figwidth 0.43*figheight 0.17*figwidth 0.04*figheight];
 posLSBAmpSet = [0.26*figwidth 0.43*figheight 0.12*figwidth 0.04*figheight];
 posSR = [0.39*figwidth 0.43*figheight 0.06*figwidth 0.04*figheight];
 posWFSel = [0.46*figwidth 0.43*figheight 0.12*figwidth 0.04*figheight];
 posBiasSet = [0.59*figwidth 0.43*figheight 0.11*figwidth 0.04*figheight];
 posselStepType = [0.71*figwidth 0.43*figheight 0.15*figwidth 0.04*figheight];
 posSuccApproxCompute = [0.87*figwidth 0.43*figheight 0.08*figwidth 0.04*figheight];

%posSetPeakValTestSig = [0.1*figwidth 0.43*figheight 0.2*figwidth 0.04*figheight];
set(hndSetPeakValTestSig,'Position',posSetPeakValTestSig);
   
%posLSBAmpSet = [0.42*figwidth 0.43*figheight 0.16*figwidth 0.04*figheight];
set(hndLSBAmpSet,'Position',posLSBAmpSet);

set(hndSRSet,'Position',posSR);

set(hndWFSel,'Position',posWFSel);

set(hndBiasSet,'Position',posBiasSet);
%posselStepType = [0.68*figwidth 0.43*figheight 0.15*figwidth 0.04*figheight];
set(selStepTypeSAR,'Position',[posselStepType(1) posselStepType(2) posselStepType(3) posselStepType(4)]);
   
%posSuccApproxCompute = [0.85*figwidth 0.43*figheight 0.08*figwidth 0.04*figheight];
set(SARCompute,'Position',posSuccApproxCompute);




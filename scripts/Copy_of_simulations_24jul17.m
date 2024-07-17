% Round of simulations started on Jul/17/2024
%
% Using dipoles at the cortical surface, normal to the cortex. EEG
% electrodes are based on 10-10 system (92 electrodes) over ICBM152.
% 
% Source patch has an effective are of 5 cm^2, never p
% laced on the sulci.
%
% Source patch has different profiles:
%     square  Unit source over a circle of radius k
%        exp  ||J_n|| = exp( - dist(n, n*) / 2*k )
%      gauss  ||J_n|| = exp( - dist(n, n*)^2 / 2*k^2 )
%       circ  ||J_n|| = sqrt( 1 - [ dist(n, n*) / k ]^2 )
%
% The reasoning for those 'profiles' is that (1) square is used on many
% papers, altought it doesn't offer smoothness in space, (2) gaussian is
% smoooth and doesn't have a heavy tail, (3) exponential is heavy-tailed,
% altough pointy, and (4) is similar to gaussian but has a clear extension.
%
% 1 single cases with SNR = 20 dB for illustrative purposes
%

%% GENERAL PARAMETERS
info = [];

% forward model
info.OGforward  = 'asa_10_10_srf_BEM';
info.OGanatomy  = 'icbm152anatomy';

info.SourceType = 'surface';

info.nTrials    = 30;
info.SNRvals    = [inf,30,20];

info.ProtocolFun   = 'Protocol04';

info.maxDepth  = 30; % unit: mm
info.maxKappa  = 10*sqrt(10/pi); % unit: mm
info.minKappa  = 10*sqrt(10/pi); % unit: mm

% for vol:  kap = 30.9 mm  ->  A = 30 cm^2
% for srf:  kap = 30.9 mm  ->  A = 30 cm^2

info.debugFigs  = false;

info.debugCent  = false;
info.debugCoord = [47.353, 18.555, 113.019];

info.print_all = false;

%% SQUARE PROFILE

info.BaseName   = 'protocol04_shape_square_vol2';
info.SourceProfile = 'square';

generator(info);
evaluator(info);
collector(info);

%% GAUSSIAN PROFILE

info.BaseName   = 'protocol04_shape_gauss_vol2';
info.SourceProfile = 'gauss';

generator(info);
evaluator(info);
collector(info);

%% EXPONENTIAL PROFILE

info.BaseName   = 'protocol04_shape_exp_vol2';
info.SourceProfile = 'exp';

generator(info);
evaluator(info);
collector(info);

%% POLYNOMIAL PROFILE

info.BaseName   = 'protocol04_shape_circ_vol2';
info.SourceProfile = 'circ';

generator(info);
evaluator(info);
collector(info);


%% RESET RESULTS
if false
  currSolver = 'MSP';
  params.(currSolver) = [];
  checklist.(currSolver) = false(nCases,1);
  checklist.tuned.(currSolver) = false;
  evaluation.(currSolver) = evaluation.(currSolver)*0;
  %
  save("params","params", '-v7.3')
  save("checklist","checklist")
  save("evaluation","evaluation")
end

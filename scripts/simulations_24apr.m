% Round of simulations 1 
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
% 500 cases over SNR (dB) levels: no-noise, 30, 20, 10, 0, -10
%

%% GENERAL PARAMETERS
info = [];

% forward model
info.OGforward  = 'asa_10_10_srf_BEM';
info.OGanatomy  = 'icbm152anatomy';

info.SourceType = 'surface';

info.nTrials    = 30;
info.SNRvals    = [Inf, 30, 20, 10, 0, -10];

info.ProtocolFun   = 'Protocol04';

info.maxDepth  = 25; % unit: mm
info.maxKappa  = 10*sqrt(5/pi); % unit: mm
info.minKappa  = 10*sqrt(5/pi); % unit: mm

% for vol:  kap = 30.9 mm  ->  A = 30 cm^2
% for srf:  kap = 30.9 mm  ->  A = 30 cm^2

info.debugFigs  = false;

info.debugCent  = false;
info.debugCoord = [47.353, 18.555, 113.019];

%% SQUARE PROFILE

info.BaseName   = 'protocol04_30_square';
info.SourceProfile = 'square';

%generator(info);
evaluator(info);
collector(info);

%% GAUSSIAN PROFILE

info.BaseName   = 'protocol04_30_gauss';
info.SourceProfile = 'gauss';

%generator(info);
evaluator(info);
collector(info);

%% EXPONENTIAL PROFILE

info.BaseName   = 'protocol04_30_exp';
info.SourceProfile = 'exp';

%generator(info);
evaluator(info);
collector(info);

%% POLYNOMIAL PROFILE

info.BaseName   = 'protocol04_30_circ';
info.SourceProfile = 'circ';

%generator(info);
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

function metric = AUROC_glo_w( meta, info, result, solution )
% Area Under Receiver-Operator Curve (AUROC), using ||J||/maxJ as label in
% the range [0,1] with the ground truth as reference.
%
% All dipoles are used for te computation.
% Using the dipole area/volume as weight, for heterogeneous dipole
% distributions.

weighted = true;
local    = false; % all points

metric = AUROC_gral( meta, result, solution, weighted, local );

end
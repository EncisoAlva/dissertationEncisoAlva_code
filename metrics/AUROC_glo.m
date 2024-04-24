function metric = AUROC_glo( meta, info, result, solution )
% Area Under Receiver-Operator Curve (AUROC), using ||J||/maxJ as label in
% the range [0,1] with the ground truth as reference.
%
% All dipoles are used for te computation.

weighted = false;
local    = false; % all points

metric = AUROC_gral( meta, result, solution, weighted, local );

end
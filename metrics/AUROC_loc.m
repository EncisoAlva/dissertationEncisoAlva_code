function metric = AUROC_loc( meta, info, result, solution )
% Area Under Receiver-Operator Curve (AUROC), using ||J||/maxJ as label in
% the range [0,1] with the ground truth as reference.
%
% Using only a neighborhood around the true center in order to alleviate
% class imbalance --too many negatives.

weighted = false;
local    = true;

metric = AUROC_gral( meta, result, solution, weighted, local );

end
function metric = AP_loc( meta, info, result, solution )
% Average Precision is the area under the Precision-Recall curve, using 
% ||J||/maxJ as label in the range [0,1] with ground truth as reference.
%
% Using only a neighborhood around the true center in order to alleviate
% class imbalance --too many negatives.

weighted = false;
local    = true;

metric = AvgPrecision_gral( meta, result, solution, weighted, local );

end
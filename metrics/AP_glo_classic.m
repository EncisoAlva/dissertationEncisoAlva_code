function metric = AP_glo_classic( meta, info, result, solution )
% Average Precision is the area under the Precision-Recall curve, using 
% ||J||/maxJ as label in the range [0,1] with ground truth as reference.
%
% All dipoles are used for te computation.

weighted = false;
local    = false; % all points

metric = AvgPrecision_classic_gral( meta, result, solution, weighted, local );

end
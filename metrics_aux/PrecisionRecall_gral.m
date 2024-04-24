function metric = AUROC_gral( meta, result, solution, weighted, local )
% Area Under Receiver-Operator Curve (AUROC), using ||J||/maxJ as label in
% the range [0,1] with the ground truth as reference.

% finding the peak
[~,IDX]    = max( solution.normJ, [], "all" );
[~, pkIDX] = ind2sub( size(solution.normJ) ,IDX);
normJpeak  = solution.normJ(:,pkIDX) / max(solution.normJ(:,pkIDX));

if local
  locDist = max( vecnorm( meta.Gridloc(result.idxDips,:) - result.TrueCent , 2, 2 ), ...
    vecnorm( meta.Gridloc(normJpeak>0.05,:) - result.TrueCent , 2, 2 ) );
  idx    = 1:meta.nDips;
  locIDX = idx( vecnorm( meta.Gridloc - result.TrueCent , 2, 2 ) < locDist );
else
  locIDX = 1:result.idxDips;
end

% score = ||J||/maxJ
% reference
Jcarrier  = zeros(1, result.nDips);
Jcarrier(result.idxDips) = result.normJshort;
score_ref = Jcarrier(locIDX);
clear Jcarrier
% estimated
score_est = normJpeak(locIDX);

% weight
if weighted
  weight = meta.GridSize(locIDX);
else
  weight = ones(1, length(locIDX));
end

% the only values of beta that actually matter
Beta  = sort(union( score_ref, score_est ), 'descend');
nBetas = length(Beta);

% TruePos, FalsePos, TrueNeg, FalseNeg
TP = zeros(1, nBetas );
FP = zeros(1, nBetas );
TN = zeros(1, nBetas );
FN = zeros(1, nBetas );
for ii = 1:nBetas
  TP(ii) = sum(weight( ( score_ref>=Beta(ii) )&&( score_est>=Beta(ii) ) ));
  FP(ii) = sum(weight( ( score_ref< Beta(ii) )&&( score_est>=Beta(ii) ) ));
  TN(ii) = sum(weight( ( score_ref< Beta(ii) )&&( score_est< Beta(ii) ) ));
  FN(ii) = sum(weight( ( score_ref>=Beta(ii) )&&( score_est< Beta(ii) ) ));
end

% ROC is parametrized by FPR vs TPR
% False Positive Rate, True Positive Rate
FPR = FP ./ ( FP + TN );
TPR = TP ./ ( TP + TN );

% area under ROC using trapezoid rule
AUROC = 0;
for ii = 2:nBetas
  AUROC = AUROC + ( FPR(ii)-FPR(ii-1) )*( TPR(ii)+TPR(ii-1) );
  % division by 2 is done at the end
end

% Area Under Receiver-Operator Curve
metric = AUROC/2;

end
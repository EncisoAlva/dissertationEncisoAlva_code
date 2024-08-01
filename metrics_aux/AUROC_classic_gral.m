function metric = AUROC_classic_gral( meta, result, solution, weighted, local )
% Area Under Receiver-Operator Curve (AUROC), using ||J||/maxJ as label in
% the range [0,1] with the ground truth as reference.

% finding the peak
[~,IDX]    = max( solution.normJ, [], "all" );
[~, pkIDX] = ind2sub( size(solution.normJ) ,IDX);
normJpeak  = solution.normJ(:,pkIDX) / max(solution.normJ(:,pkIDX));

% local: use only sources with magnitude higher than 10% of maxJ
hard_max = 10*sqrt(30/pi);
if local
  locDist = min( [...
    max(vecnorm( meta.Gridloc(result.data.idxShort,:) - result.data.TrueCent , 2, 2 )), ...
    max(vecnorm( meta.Gridloc(normJpeak>0.1,:) - result.data.TrueCent , 2, 2 ) ), ...
    hard_max] );
  idx    = 1:meta.nGridDips;
  switch meta.info.SourceType
    case 'volume'
      locIDX = idx( vecnorm( meta.Gridloc - result.data.TrueCent , 2, 2 ) <= locDist );
    case 'surface'
      [~,GraphDist] = shortestpathtree(meta.asGraph, result.idxCent, idx );
      locIDX = idx( GraphDist <= locDist );
  end
else
  locIDX = 1:1:meta.nGridDips;
end

% score = ||J||/maxJ
% reference
Jcarrier  = zeros(1, meta.nGridDips);
Jcarrier(result.data.idxShort) = result.data.normJshort;
score_ref0 = Jcarrier(locIDX);
score_ref = max(score_ref0(:))*( score_ref0 > max(score_ref0(:))*.5 );
clear Jcarrier
% estimated
score_est = normJpeak(locIDX)';

% weight
if weighted
  weight = meta.GridVolume(locIDX);
else
  weight = ones(1, length(locIDX));
end

% the only values of beta that actually matter
Beta  = sort(unique( union( score_ref, score_est )), 'descend');
nBetas = length(Beta);

% TruePos, FalsePos, TrueNeg, FalseNeg
TP = zeros(1, nBetas );
FP = zeros(1, nBetas );
TN = zeros(1, nBetas );
FN = zeros(1, nBetas );
for ii = 1:nBetas
  TP(ii) = sum(weight( and( score_ref>=Beta(ii), score_est>=Beta(ii) ) ));
  FP(ii) = sum(weight( and( score_ref< Beta(ii), score_est>=Beta(ii) ) ));
  TN(ii) = sum(weight( and( score_ref< Beta(ii), score_est< Beta(ii) ) ));
  FN(ii) = sum(weight( and( score_ref>=Beta(ii), score_est< Beta(ii) ) ));
end

% debug figure
if meta.info.debugFigs
  if local
    lab_loc = 'Local';
  else
    lab_loc = 'Non-local';
  end
  if weighted
    lab_we = 'Weighted';
  else
    lab_we = 'Non-weighted';
  end
  figure()
  plot(Beta, TP, Beta, FP, Beta, TN, Beta, FN)
  xlabel('Classification threshold, beta')
  legend('True Positives', 'False Positives', 'True Negatives', 'False Negatives')
  if weighted
    switch meta.info.SourceType
      case 'surface'
        ylabel('Area [mm^2]')
      case 'volume'
        ylabel('Voluma [mm^2]')
    end
  else
    ylabel('Dipole count')
  end
  subtitle([lab_loc, ', ', lab_we])
end

% ROC is parametrized by FPR vs TPR
% False Positive Rate, True Positive Rates
FPR = FP ./ ( FP + TN +1);
TPR = TP ./ ( TP + TN +1);

FPR(1)   = 0;
FPR(end) = 1;
TPR(1)   = 0;
TPR(end) = 1;

if meta.info.debugFigs
  figure()
  tiledlayout(1,2)
  %
  nexttile
  plot(Beta, TPR, Beta, FPR)
  xlabel('Classification threshold')
  legend('True Positive Rate', 'False Positive Rate')
  ylabel('Unitless')
  %
  nexttile
  plot(FPR, TPR)
  xlabel('False Positive Rate')
  ylabel('True Positive Rate')
  subtitle([lab_loc, ', ', lab_we])
end

if meta.info.debugFigs
  figure()
  plot(FPR, TPR)
  grid on
  xlabel('FPR = FP/(FP+TN)', 'Interpreter','latex')
  ylabel('TPR = TP/(TP+FN)', 'Interpreter','latex')
  %
  hold on
  for bb = 0:0.2:1
    [~, ii] = min(abs(bb-Beta));
    text(FPR(ii), TPR(ii),  ...
      ['$\beta = {',num2str(bb),'}$'], 'Interpreter','latex')
  end
  for bb = 0:.05:1
    [~, ii] = min(abs(bb-Beta));
    scatter(FPR(ii), TPR(ii), 15,'blue','filled')
  end
  title('Receiver-Operator Curve')
  set(gcf,'color','w');
  set(gca,'LooseInset',get(gca,'TightInset'))
  fig = gcf;
  fig.Units = 'inches';
  fig.OuterPosition = [0 0 3 3]*1.5;
  exportgraphics(gcf,'ROC.pdf','Resolution',600)
end

% area under ROC using trapezoid rule
AUROC = 0;
for ii = 2:nBetas
  if ~isnan(FPR(ii)-FPR(ii-1)) && ~isnan(TPR(ii)-TPR(ii-1))
    AUROC = AUROC + ( FPR(ii)-FPR(ii-1) )*( TPR(ii)+TPR(ii-1) );
    % division by 2 is done at the end
  end
end

% Area Under Receiver-Operator Curve
metric = AUROC/2;

end
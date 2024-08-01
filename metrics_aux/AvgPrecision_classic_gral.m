function metric = AvgPrecision_classic_gral( meta, result, solution, weighted, local )
% Average Precision is the area under the Precision-Recall curve, using 
% ||J||/maxJ as label in the range [0,1] with ground truth as reference.

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
  xlabel('Classification threshold')
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

% Average Precision is the area under the Precision-Recall curve
rec = TP ./ ( TP + FN +1 ); % the +1 is to avoid NaNs
pre = TP ./ ( TP + FP +1 );

pre(end) = pre(end-1);

if meta.info.debugFigs
  figure()
  tiledlayout(1,2)
  %
  nexttile
  plot(Beta, pre, Beta, rec)
  xlabel('Classification threshold')
  legend('Precision', 'Recall')
  ylabel('Unitless')
  %
  nexttile
  plot(rec, pre)
  ylabel('Precision')
  xlabel('Recall')
  subtitle([lab_loc, ', ', lab_we])
end

if meta.info.debugFigs
  figure()
  plot(rec, pre)
  grid on
  xlabel('Recall = TP/(TP+FN)', 'Interpreter','latex')
  ylabel('Precision = TP/(TP+FP)', 'Interpreter','latex')
  xlim([0,1])
  ylim([0,1])
  %
  hold on
  for bb = 0:0.2:1
    [~, ii] = min(abs(bb-Beta));
    text(rec(ii), pre(ii),  ...
      ['$\beta = {',num2str(bb),'}$'], 'Interpreter','latex')
  end
  for bb = 0:.05:1
    [~, ii] = min(abs(bb-Beta));
    scatter(rec(ii), pre(ii), 15,'blue','filled')
  end
  title('Precision-Recall Curve')
  set(gcf,'color','w');
  set(gca,'LooseInset',get(gca,'TightInset'))
  fig = gcf;
  fig.Units = 'inches';
  fig.OuterPosition = [0 0 3 3]*1.5;
  exportgraphics(gcf,'PRC.pdf','Resolution',600)
end

% area under ROC using trapezoid rule
AP = 0;
for ii = 2:nBetas
  if ~isnan(rec(ii)-rec(ii-1)) && ~isnan(pre(ii)-pre(ii-1))
    AP = AP + ( rec(ii)-rec(ii-1) )*( pre(ii)+pre(ii-1) );
    % division by 2 is done at the end
  end
end

% Average Precision
metric = AP/2;

end
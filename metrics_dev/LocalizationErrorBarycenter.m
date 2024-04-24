function metric = LocalizationErrorWeighted( meta, result, solution )
% Source location is defined as the a weighted sum of dipole locations,
% with the weight being their magnitude
%    (src loc) = (sum_n loc_n * ||J_n||_2) / (sum_n ||J_n||_2)

% time averaging
Javg = movmean( solution.J, 3 );
switch meta.type
  case 'surface'
    % nothing else
  case 'volume'
    n    = size(Javg,1)/3;
    tmp  = Javg;
    Javg = zeros(n, size(Javg,2));
    for i = 1:n
      Javg(i,:) = vecnorm(tmp(3*(i-1)+(1:3), :), 2, 1);
    end
end

% identify time with the maximum peak
[~,IDX]    = max( Javg, [], "all" );
[peakIDX, ~] = ind2sub( size(solution.J) ,IDX);

metric = norm( Javg(peakIDX,:)*meta.Gridloc(peakIDX,:)/sum(Javg(peakIDX,:)) - result.data.TrueBarycent )*1000;

%metric = norm( meta.Gridloc(peakIDX,:) - result.data.TrueBarycent )*1000;

end
function metric = SpaDis_gral( meta, info, result, solution, p )
% Spatial Dispersion of order p, defined as
%  SD_p = ( sum_n [ dist(n,n*) ||J_n||_2 ]^p / sum_n [ ||J_n||_2 ]^p )^1/p
%
% The case p=2 is the one proposed by Molins et al (2008), but I think that
% p=1 makes more sense as a weighted sum.
% Also, I think that using the estimated barycenter instead of the true
% barycenter would be useful for real data.

% finding the peak
[~,IDX]    = max( solution.normJ, [], "all" );
[~, pkIDX] = ind2sub( size(solution.normJ) ,IDX);
normJpeak  = solution.normJ(:,pkIDX);

% distance to true barycenter
dist_true = vecnorm( meta.Gridloc - result.data.TrueCent, 2, 2 );

% Spatial Dispersion of order p
metric = ( sum( (dist_true.*normJpeak).^p ) / sum( normJpeak.^p ) )^(1/p);

end
function metric = DLE1( meta, info, result, solution )
% Source location is defined as the a weighted sum of dipole locations,
% with the weight being their magnitude
%    (src loc) = (sum_n loc_n * ||J_n||_2) / (sum_n ||J_n||_2)

% finding the peak
[~,IDX]    = max( solution.normJ, [], "all" );
[~, pkIDX] = ind2sub( size(solution.normJ) ,IDX);

% center of mass of J at the peak
estCenter = solution.normJ(:,pkIDX)' * meta.Gridloc / sum( solution.normJ(:) );

% reference center of mass
refCenter = result.data.TrueCent;

% Distance Localisation Error for one single point
metric = norm( estCenter - refCenter );

end
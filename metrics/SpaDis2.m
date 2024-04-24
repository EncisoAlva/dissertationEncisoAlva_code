function metric = SpaDis2( meta, info, result, solution )
% Spatial Dispersion of order 2, defined as
%  SD_1 = ( sum_n [ dist(n,n*) ||J_n||_2 ]^2 / sum_n [ ||J_n||_2 ]^2 )^1/2
%
% The case p=2 is the one proposed by Molins et al (2008), but I think that
% p=1 makes more sense as a weighted sum.

% Particular case
metric = SpaDis_gral( meta, info, result, solution, 2 );

end
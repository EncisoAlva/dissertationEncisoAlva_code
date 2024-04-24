function metric = HalfMax( meta, info, result, solution )
% Tag the dipoles such that J_n > 0.5 * maxJ. Report the surface/volume of
% the region defined by those dipoles.

% finding the peak
[~,IDX]    = max( solution.normJ, [], "all" );
[~, pkIDX] = ind2sub( size(solution.normJ) ,IDX);
Jest_pk    = solution.normJ(:,pkIDX);

% Half-Max Size
metric = sum( meta.GridVolume( Jest_pk > 0.5*max(Jest_pk) ) );

end
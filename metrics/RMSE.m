function metric = RMSE( meta, info, result, solution )
% Relative Mean Square Error
%    || J^ - J ||_2 / || J ||_2

% finding the peak
[~,IDX]    = max( solution.normJ, [], "all" );
[~, pkIDX] = ind2sub( size(solution.normJ) ,IDX);

% ground truth
estJ = solution.J(:,pkIDX);
switch info.SourceType
  case 'surface'
    refJ = zeros( meta.nGridDips, 1 );
    refJ( result.data.idxShort ) = result.data.Jshort;
    %
    relRMSE = norm( estJ - refJ, 2 ) / norm( refJ, 2 );
  case 'volume'
    refJ = zeros( meta.nDips*3, 1 );
    refJ( (result.idxShort*3)+[1,2,3] ) = result.Jshort;
    %
    relRMSE = dip_norm( estJ - refJ, 2 ) / norm( refJ, 2 );
end

% Distance Localisation Error for one single point
metric = relRMSE;

end
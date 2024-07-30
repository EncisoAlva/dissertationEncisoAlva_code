function metric = funMSE( meta, info, result, solution )
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
  case 'volume'
    refJ = zeros( meta.nGridDips*3, 1 );
    refJ( result.data.idxShortG ) = result.data.Jshort;
end
%
funMSE = norm( meta.Leadfield*(estJ - refJ), 2)  / norm( meta.Leadfield*refJ, 2 );

% Distance Localisation Error for one single point
metric = funMSE;

end
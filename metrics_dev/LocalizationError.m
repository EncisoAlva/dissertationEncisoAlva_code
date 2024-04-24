function metric = LocalizationError( meta, result, solution )
% Source location is defined as the dipole with max magnitude,
%    (src loc) = argmax_n ||J_n||_2

% time averaging
Javg = movmean( solution.J, 3 );
switch meta.type
  case 'surface'
    Javg = abs(Javg);
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

metric = norm( meta.Gridloc(peakIDX,:) - result.data.TrueBarycent )*1000;

end
function metric = HalfMax( meta, result, solution )

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
[Jmax,IDX]    = max( Javg, [], "all" );
[~,peakTIME] = ind2sub( size(solution.J) ,IDX);

metric = sum(meta.GridVolume( (Javg(:,peakTIME)/Jmax) > 0.5 ))*(1000^2);

end
function solution = SingleRegionPrior( meta, info, result, pars )
% TODO add description (optional)
%
solution = [];

% solution per se
parTic = tic;

% hot start: the MNE solution
oldJ = meta.Leadfield' * ...
  pinv( meta.Leadfield * meta.Leadfield' + pars.lambda* eye(pars.M) ) * result.data.Y;
switch info.SourceType
  case 'surface'
    oldJnorm = abs(oldJ).^2;
  case 'volume'
    oldJnorm = dip_norm(oldJ).^2;
end

%initial partition of R0
g0 = pars.g0;
g1 = pars.g1;
N0 = pars.N0;

% loop
counter = 0; err = Inf;
while (counter < 20) && ( err>10^-5 )
  counter = counter + 1;
  for ii = 1:5
    % identification of regions based on a threshold
    thr = ( N0*g0 + (pars.N-N0)*g1 ) / pars.N;
    if (thr < min(oldJnorm)) || (thr > max(oldJnorm))
      thr = ( min(oldJnorm) + max(oldJnorm) )/2;
    end
    R0  = ( oldJnorm <= thr );
    R1  = ( oldJnorm >  thr );
    N0  = sum(R0);
    % copmutation of parameters from each region
    g0  = mean( oldJnorm(R0) );
    g1  = mean( oldJnorm(R1) );
  end
  if info.debugFigs
    figure()
    tiledlayout(2,2)
    nexttile
    trisurf(meta.Cortex.Faces, ...
      meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
    hold on
    scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
      40, oldJnorm*120,'filled')
    colormap("parula")
    scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
      200, 'red','filled')
    title("Magnitude of estimated sources")
    b = colorbar;
    b.Label.String = 'Unitless; range=[0,120]';
    %
    nexttile
    trisurf(meta.Cortex.Faces, ...
      meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
    hold on
    scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
      40, mean(oldJnorm)*R1*120,'filled')
    colormap("parula")
    scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
      200, 'red','filled')
    title("Active region")
    colorbar;
    %
    nexttile
    pd = makedist('Gamma', 'a', 1/2, 'b', 2); % chi2 1dof
    qqplot(oldJnorm(R0)/g0, pd);
    %
    nexttile
    qqplot(oldJnorm(R1)/g1, pd);
  end

  % averaging operator
  switch info.SourceType
    case 'surface'
      Lk = zeros(pars.N, 1);
      Lk( R1 ) = 1;
      nu = pars.N;
    case 'volume'
      Lk = zeros(pars.N*3, 1);
      Lk( (kron(R1, [1,1,1]'))>0 ) = 1;
      nu = pars.N*3;
  end
  Hs = sparse( sparse(g0*eye(nu)) + g1*(1/sum(R1))*sparse((Lk*Lk')) );

  % new solution
  newJ = Hs * meta.Leadfield' * ...
    pinv( meta.Leadfield*Hs*meta.Leadfield' + pars.lambda*eye(pars.M) ) * result.data.Y;
  err = norm( newJ-oldJ, 'fro' );

  % update for next iteration, if any
  oldJ = newJ;
  switch info.SourceType
    case 'surface'
      oldJnorm = abs(oldJ).^2;
    case 'volume'
      oldJnorm = dip_norm(oldJ).^2;
  end
  maxJ = max(oldJnorm);

  % print partial results
  fprintf("Iter %2d. Err = %3.3d. (0.00): %5d dips. (0.05) : %4d dips.  (0.50) : %3d dips.\n", ...
    counter, err, sum(oldJnorm~=0), sum(oldJnorm>0.05*maxJ), sum(oldJnorm>0.5*maxJ))
end

% solution per se
solution.J = oldJ;
solution.algTime = toc(parTic);

% norm of J
switch meta.Type
  case 'surface'
    solution.normJ = abs( solution.J );
  case 'volume'
    solution.normJ = dip_norm( solution.J );
end

% stop timer
solution.algTime = toc(parTic);

end
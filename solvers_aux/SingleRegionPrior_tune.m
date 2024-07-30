function pars = SingleRegionPrior_tune( meta, info, result )
% TODO add description (optional)
%

% start timer
parTic = tic;
pars   = [];

% general values
pars.M = meta.nChans;
pars.N = meta.nGridDips;
pars.T = size(result.data.time, 2);

% hot start: the MNE solution
oldJ = meta.Leadfield' * ...
  pinv( meta.Leadfield * meta.Leadfield' + median(meta.S)^2* eye(pars.M) ) * result.data.Y;
switch info.SourceType
  case 'surface'
    oldJnorm = abs(oldJ).^2;
  case 'volume'
    oldJnorm = dip_norm(oldJ).^2;
end

%initial partition of R0
g0 = min( oldJnorm );
g1 = max( oldJnorm );
N0 = max( sum( oldJnorm < (g0+g1)/2 ), pars.N/2 );
best_lambda = median(meta.S)^2*( g1 );

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
    %
    if info.debugFigs
      if ii == 1
        figure()
        histogram(oldJnorm)
        hold on
      else
        xline(thr, '-', {num2str(ii),' iter'})
      end
    end
    % copmutation of parameters from each region
    g0  = mean( oldJnorm(R0) );
    g1  = mean( oldJnorm(R1) );
  end
  if info.debugFigs
    figure()
    tiledlayout(2,2)
    nexttile([1,2])
    histogram(oldJnorm)
    xline(thr, '-', {'Threshold'})
    %
    nexttile
    histogram(oldJnorm(R0))
    nexttile
    histogram(oldJnorm(R1))
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
      40, R1*120,'filled')
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
  Hs = sparse(g0*eye(nu)) + sparse(g1*(1/sum(R1))*(Lk*Lk'));
  %As = sparse( sparse(eye(pars.N)) - sparse((1/sum(R1))*(Lk*Lk')) );

  % hyperparameter tuning via Generalized Cross-Validation
  % provisional, just to obtain the region properly
  scale  = 10;
  for iter = 1:3
    % try many values for lambda, compute GCV value for each, get the min
    lambdas = best_lambda * (2.^( (-scale):(scale/10):scale ));
    Gs     = zeros( size(lambdas) );
    for q = 1:length(lambdas)
      lambda = lambdas(q);
      %Gs(q) = SingleRegionPrior_Lcurve( meta, info, result, pars, lambda, Hs, g0, g1, As );
      Gs(q) = SingleRegionPrior_GCV( meta, info, result, pars, lambda, Hs );
      %fprintf("Lambda = %2.3d. L-coord = %2.3d.\n", lambda, Gs(q))
    end
    [~, idx]   = min(Gs);
    best_lambda = lambdas(idx);
    %
    % if not on the border, reduce scale; else, increase it
    if (1<idx) && (idx<length(Gs))
      scale = scale/10;
    else
      scale = scale*10;
    end
  end
  best_lambda = max(best_lambda, 0.001);

  % new solution
  newJ = Hs * meta.Leadfield' * ...
    pinv( meta.Leadfield*Hs*meta.Leadfield' + best_lambda*eye(pars.M) ) * result.data.Y;
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

% X ~ chi(1)
%   E[   X ] = 1
% Var(   X ) = 2
%   E[ k*X ] = k
% Var( k*X ) = 2*k^2

% hyperparameter tuning via Generalized Cross-Validation
% actual tuning for lambda
% re-initialize to escape local minima
best_lambda = median(meta.S)^2;
scale  = 10;
for iter = 1:6
  % try many values for lambda, compute GCV value for each, get the min
  lambdas = best_lambda * (2.^( (-scale):(scale/10):scale ));
  Gs     = zeros( size(lambdas) );
  for q = 1:length(lambdas)
    lambda = lambdas(q);
    %Gs(q) = SingleRegionPrior_Lcurve( meta, info, result, pars, lambda, Hs, g0, g1, As );
    Gs(q) = SingleRegionPrior_GCV( meta, info, result, pars, lambda, Hs );
  end
  [~, idx]   = min(Gs);
  best_lambda = lambdas(idx);
  %
  % if not on the border, reduce scale; else, increase it
  if (1<idx) && (idx<length(Gs))
    scale = scale/10;
  else
    scale = scale*10;
  end
end
pars.lambda = max(best_lambda, 0.001);
pars.g0     = g0;
pars.g1     = g1;
pars.N0     = N0;

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for SingleRegionPrior solver.\n Optimal lambda: ")
disp(pars.lambda)
fprintf("\n")

end
function pars = SingleRegionPriorBay_tune( meta, info, result )
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
J = meta.Leadfield' * ...
  pinv( meta.Leadfield * meta.Leadfield' + median(meta.S)^2* eye(pars.M) ) * result.data.Y;
switch info.SourceType
  case 'surface'
    Jnorm = abs(J).^2;
  case 'volume'
    Jnorm = dip_norm(J).^2;
end

%initial partition of R0
prob_R0 = ones( pars.N, 1 )*0.5;
prob_R1 = ones( pars.N, 1 )*0.5;
prob_R0( Jnorm <  mean(Jnorm) ) = .99;
prob_R1( Jnorm >  mean(Jnorm) ) = .99;
%
if ( sum( (prob_R0> prob_R1)*1 ) < 10 )||( sum( (prob_R0<=prob_R1)*1 ) < 10 )
  %break
  ab0_ = [1,1];
  ab1_ = [1,1];
end
ab0_ = gamfit( Jnorm(prob_R0> prob_R1) );
ab1_ = gamfit( Jnorm(prob_R0<=prob_R1) );
if any(isnan(ab0_)) || any(isnan(ab1_))
  %break
  ab0 = [1,1];
  ab1 = [1,1];
else
  ab0 = ab0_;
  ab1 = ab1_;
end
%
best_lambda = median(meta.S)^2;

% loop
counter = 0; err = Inf;
while (counter < 20) && ( err>10^-5 )
  counter = counter + 1;
  % finding an appropriate partition of R0 and R1
  for ii = 1:100
    poster = gampdf(Jnorm, ab0(1),ab0(2)).*prob_R0 + gampdf(Jnorm, ab1(1),ab1(2)).*prob_R1;
    prob_R0_new = ( gampdf(Jnorm, ab0(1),ab0(2)).*prob_R0 ) ./ poster;
    prob_R1_new = ( gampdf(Jnorm, ab1(1),ab1(2)).*prob_R1 ) ./ poster;
    if info.debugFigs
      figure()
      tiledlayout(2,2)
      nexttile
      histogram(prob_R0)
      title('Prob(n in R_0)^{(k)}')
      nexttile
      histogram(prob_R0_new)
      title('Prob(n in R_0)^{(k+1)}')
      nexttile
      histogram(prob_R1)
      title('Prob(n in R_1)^{(k)}')
      nexttile
      histogram(prob_R1_new)
      title('Prob(n in R_1)^{(k+1)}')
      %
      figure()
      semilogx( Jnorm, prob_R0_new, '.' )
      hold on
      semilogx( Jnorm, prob_R1_new, '.' )
      xlabel('|| J ||_2^2')
      %
      figure()
      tiledlayout(2,2)
      nexttile
      histogram(Jnorm(prob_R0_new> prob_R1_new))
      nexttile
      histogram(Jnorm(prob_R0_new<=prob_R1_new))
      nexttile
      pd0 = makedist('Gamma', 'a', ab0(1), 'b', ab0(2));
      qqplot(Jnorm(prob_R0_new> prob_R1_new), pd0);
      nexttile
      pd1 = makedist('Gamma', 'a', ab1(1), 'b', ab1(2));
      qqplot(Jnorm(prob_R0_new<=prob_R1_new), pd1);
    end
    prob_R0 = prob_R0_new;
    prob_R1 = prob_R1_new;
    %
    if ( sum( (prob_R0> prob_R1)*1 ) < 10 )||( sum( (prob_R0<=prob_R1)*1 ) < 10 )
      break
    end
    ab0_ = gamfit( Jnorm(prob_R0> prob_R1) );
    ab1_ = gamfit( Jnorm(prob_R0<=prob_R1) );
    if any(isnan(ab0_)) || any(isnan(ab1_))
      break
    else
      ab0 = ab0_;
      ab1 = ab1_;
    end
  end
  disp([ab0, ab1])

  % averaging operator
  R1 = (prob_R1>prob_R0);
  %R1 = (Jnorm > gaminv( .5, ab1(1),ab1(2)) );
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
  Hs = sparse( ab0(1)*(ab0(2)^2)*eye(nu)) + ...
      sparse(ab1(1)*(ab1(2)^2)*(1/sum(R1))*(Lk*Lk'));
  %As = sparse( sparse(eye(nu)) - sparse((1/sum(R1))*(Lk*Lk')) );

  % hyperparameter tuning via Generalized Cross-Validation
  % provisional, just to obtain the region properly
  scale  = 10;
  for iter = 1:3
    % try many values for lambda, compute GCV value for each, get the min
    lambdas = best_lambda * (2.^( (-scale):(scale/10):scale ));
    Gs     = zeros( size(lambdas) );
    for q = 1:length(lambdas)
      lambda = lambdas(q);
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
  err = norm( newJ-J, 'fro' );

  % update for next iteration, if any
  J = newJ;
  switch info.SourceType
    case 'surface'
      Jnorm = abs(J).^2;
    case 'volume'
      Jnorm = dip_norm(J).^2;
  end
  maxJ = max(Jnorm);

  % print partial results
  fprintf("Iter %2d. Err = %3.3d. (0.00): %5d dips. (0.05) : %4d dips.  (0.50) : %3d dips.\n", ...
    counter, err, sum(Jnorm~=0), sum(Jnorm>0.05*maxJ), sum(Jnorm>0.5*maxJ))
end

% some debug figures for visualization
if info.debugFigs
  figure()
  tiledlayout(2,3)
  nexttile
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
    40, Jnorm/max(Jnorm),'filled')
  colormap("parula")
  clim([0,1])
  scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
    200, 'red','filled')
  title("Normalized magnitude of estimated sources")
  b = colorbar;
  b.Label.String = 'Unitless on [0,1]';
  %
  nexttile
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
    40, prob_R0,'filled')
  colormap("parula")
  clim([0,1])
  scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
    200, 'red','filled')
  title("Probabilities for R0")
  colorbar;
  %
  nexttile
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
    40, prob_R1,'filled')
  colormap("parula")
  clim([0,1])
  scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
    200, 'red','filled')
  title("Probabilities for R1")
  colorbar;
  %
  nexttile
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  %scatter3(meta.Gridloc(R1,1), meta.Gridloc(R1,2), meta.Gridloc(R1,3), ...
  %  40, gampdf(Jnorm(R1), ab1(1),ab1(2)),'filled')
  scatter3(meta.Gridloc(R1,1), meta.Gridloc(R1,2), meta.Gridloc(R1,3), ...
    40, gamcdf(Jnorm(R1), ab1(1),ab1(2)),'filled')
  colormap("parula")
  clim([0,1])
  scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
    200, 'red','filled')
  title("pdf_{gam}( ||J|| )")
  colorbar;
  %
  nexttile
  pd0 = makedist('Gamma', 'a', ab0(1), 'b', ab0(2));
  qqplot(Jnorm(prob_R0_new> prob_R1_new), pd0);
  nexttile
  pd1 = makedist('Gamma', 'a', ab1(1), 'b', ab1(2));
  qqplot(Jnorm(prob_R0_new<=prob_R1_new), pd1);
end

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

% save for using with others
pars.lambda = max(best_lambda, 10^-5);
pars.ab0    = ab0;
pars.ab1    = ab1;

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for Single Region Prior (bayesian) solver.\n Optimal lambda: ")
disp(pars.lambda)
fprintf("\n")

end

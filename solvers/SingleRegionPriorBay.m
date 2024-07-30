function solution = SingleRegionPriorBay( meta, info, result, pars )
% TODO add description (optional)
%
solution = [];

% solution per se
parTic = tic;

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
ab0 = pars.ab0;
ab1 = pars.ab1;
%
best_lambda = pars.lambda;

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
  %As = sparse( sparse(eye(pars.N)) - sparse((1/sum(R1))*(Lk*Lk')) );

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

% at this point, the solution was already computed
solution.algTime = toc(parTic);

%solution.J = diag(meta.ColumnNorm.^-1) * J;
solution.J = J;

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
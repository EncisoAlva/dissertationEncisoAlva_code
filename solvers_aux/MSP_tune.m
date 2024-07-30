function pars = MSP_tune( meta, info, result )
% TODO add description (optional)
%

% start timer
parTic = tic;
pars   = [];

% general values
pars.m  = meta.nChans;
pars.n  = meta.nGridDips;
pars.t  = size(result.data.time, 2);

% find an appropriate number of eigenvalues
[pars.Unorm,pars.Snorm,~] = svd(meta.LeadfieldColNorm,'vector');
tmp = cumsum(pars.Snorm) / sum(pars.Snorm);
pars.s = find( tmp > 0.9, 1 );

% project measurements over the first s eigenvectors of G
Us = pars.Unorm(:, (1:pars.s));
Ys = Us * Us' * result.data.Y / norm( result.data.Y );

% project G into the projection of Y
Ps = Ys * pinv( Ys'*Ys ) * Ys';

% Activation Probability Map (APM) is defined as follows
D = zeros( meta.nGridDips, 1 );
switch info.SourceType
  case 'surface'
    for ii = 1:meta.nGridDips
      D(ii) = norm( Ps * meta.LeadfieldColNorm(:,ii), 2 )^2;
    end
  case 'volume'
    for ii = 1:meta.nGridDips
      D(ii) = 0;
      for tt = 1:3
          D(ii) = D(ii) + norm( Ps * meta.LeadfieldColNorm(:,3*(ii-1)+tt), 2 )^2/3;
      end
    end
end
pars.APM = D;
if info.debugFigs
  figure()
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
    40, D*120,'filled')
  colormap("parula")
   scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
    200, 'red','filled')
  title("Activation Probability Map")
  b = colorbar;
  b.Label.String = 'Unitless; range=[0,120]';
end

% fitting a gamma distribution
% D ~ Gamma(a,b)
%pars.gamma_b = var(D) / mean(D);
%pars.gamma_a = mean(D) / pars.gamma_b;
[mle,~]  = gamfit(D);
pars.gamma_a = mle(1);
pars.gamma_b = mle(2);

% debug figures
if info.debugFigs
  figure()
  histogram(D, 'Normalization', 'pdf')
  hold on
  plot((0:0.01:1), gampdf( (0:0.01:1), pars.gamma_a, pars.gamma_b))
  title('Fitting APM as a gamma distribution')
  xlabel('Activation Probability Map (APM) value')
  ylabel('Histogram, normalized to compare with pdf')
end

% APM=1 means that a dipole is 'active', with D ~ gamma(a,b)
% I will decide a dipoles is active if P(APM>0)>threshold
pars.APMthreshold_non0 = icdf('gamma', 0.05, pars.gamma_a, pars.gamma_b );
%
pars.APMthreshold_top1 = icdf('gamma', 0.90, pars.gamma_a, pars.gamma_b );
idx = 1:meta.nGridDips;
pars.ActiveDips = idx( pars.APM > pars.APMthreshold_non0 );

% another thing to do with APM is to create weights
pars.W   = 1 - pars.APM *.99; % robustness(?)
switch info.SourceType
  case 'surface'
    % nothing else
  case 'volume'
    pars.W = kron( pars.W, [1,1,1]' );
end
pars.GWG = meta.LeadfieldColNorm * diag( pars.W.^(-2) ) * meta.Leadfield';

% hyperparameter tuning via Generalized Cross-Validation
% starting at median eigenvalue
best_alpha = median(pars.Snorm)^2;
scale  = 10;
for iter = 1:6
  % try many values for alpha, compute GCV value for each, get the min
  alphas = best_alpha * (2.^( (-scale):(scale/10):scale ));
  Gs     = zeros( size(alphas) );
  for q = 1:length(alphas)
    alpha = alphas(q);
    Gs(q) = MSP_GCV( meta, info, result, pars, alpha );
  end
  [~, idx]   = min(Gs);
  best_alpha = alphas(idx);
  %
  % if not on the border, reduce scale; else, increase it
  if (1<idx) && (idx<length(Gs))
    scale = scale/10;
  else
    scale = scale*10;
  end
end
pars.alpha  = max(best_alpha, 0.001);
pars.kernel = diag( pars.W.^(-2) ) * meta.LeadfieldColNorm' * pinv( pars.GWG + pars.alpha*eye(pars.m) );

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for MSP solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
function pars = RegionPriors_tune( meta, info, result )
% TODO add description (optional)
%

% start timer
parTic = tic;
pars   = [];

% general values
pars.M = size(meta.Leadfield, 1);
pars.N = size(meta.Leadfield, 2);
pars.T = size(result.data.time, 2);

% hot start: the MNE solution
Jold = meta.Leadfield' * ...
  pinv( meta.Leadfield * meta.Leadfield' + median(meta.S)^2* eye(pars.m) ) * result.data.Y;

% loop
iter = 0;
while iter < 100
  %
end





pars.GWG = meta.Leadfield * diag( pars.W.^(-2) ) * meta.Leadfield';

% hyperparameter tuning via Generalized Cross-Validation
% starting at median eigenvalue
best_alpha = median(meta.S)^2;
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
pars.kernel = diag( pars.W.^(-2) ) * meta.Leadfield' * pinv( pars.GWG + pars.alpha*eye(pars.m) );

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for MSP solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
function pars = sLORETA_tune( meta, info, result )
% TODO add description (optional)
%

% start timer
parTic = tic;
pars   = [];

% general values
pars.M    = size(meta.Leadfield, 1);
pars.N    = size(meta.Leadfield, 2);
pars.T    = size(result.data.time, 2);
pars.H    = eye(pars.M) - (ones(pars.M)/pars.M); % averaging operator
pars.HGGH = pars.H * meta.Leadfield * meta.Leadfield' * pars.H;

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
    Gs(q) = sLORETA_GCV( meta, info, result, pars, alpha );
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
pars.kernel = meta.Leadfield' * pars.H * pinv( pars.HGGH + pars.alpha*eye(pars.M) );

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for sLORETA solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
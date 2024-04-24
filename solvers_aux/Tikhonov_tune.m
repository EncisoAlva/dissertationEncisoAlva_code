function pars = Tikhonov_tune( meta, info, result )
% Tuning via Generalized Cross-Validation for wMNE, and additional
% initialization routines.
%
% Weighten Minimum-Norm Estimator (wMNE), follows the basic Tikonov
% regularized estimation
%   J^ = argmin_J || G*J-Y ||^2_F + alpha || W*J ||^2_F
% with W the weight induced by column-normalization of G and ||*||_F is the
% Frobenius norm.
%

% start timer
parTic = tic;
pars   = [];

% general values
pars.m = size(meta.Leadfield, 1);
pars.n = size(meta.Leadfield, 2);
pars.r = min(pars.m, pars.n);
pars.t = size(result.data.time, 2);

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
    Gs(q) = Tikhonov_GCV( meta, result, pars, alpha );
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
pars.alpha = max(best_alpha, 0.001);

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for wMNE solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
function pars = wMNE_tune( meta, info, result )
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
pars.m = size(meta.LeadfieldColNorm, 1);
pars.n = size(meta.LeadfieldColNorm, 2);
pars.r = min(pars.m, pars.n);
pars.t = size(result.data.time, 2);

%% hyperparameter tuning via Generalized Cross-Validation
% starting at median eigenvalue

alphas = 10.^(-5:0.1:5);
Gs     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Gs(q) = wMNE_GCV( meta, result, pars, alpha );
end

[~, idx]   = min(Gs);
best_alpha = alphas(idx);
GCV_alpha  = best_alpha;

pars.alpha = max(GCV_alpha, 0.001);

pars.kernel = meta.LeadfieldColNorm' * pinv( eye(pars.m) * pars.alpha + meta.LeadfieldColNorm * meta.LeadfieldColNorm' );

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for wMNE solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
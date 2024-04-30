function solution = wMNE( meta, info, result, pars )
% Weighten Minimum-Norm Estimator (wMNE), follows the basic Tikonov
% regularized estimation
%   J^ = argmin_J || G*J-Y ||^2_F + alpha || W*J ||^2_F
% with W the weight induced by column-normalization of G and ||*||_F is the
% Frobenius norm.
%
% The same is achieved by column-normalizing G and unit weight for J.
%
% The parameter alpha is found using Generalized Cross Validation.
%
solution = [];

% initialize timer
parTic   = tic;

% kernel
K = meta.Leadfield' * pinv( eye(pars.m) * pars.alpha + meta.Leadfield * meta.Leadfield' );

% solution
solution.J = K * result.data.Y;

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
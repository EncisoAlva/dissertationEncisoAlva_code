function solution = Tikhonov( meta, info, result, pars )
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

% intialize
J = zeros( pars.n, pars.t );

% kernel is not computed, only solution using SVD decompostoin
parTic   = tic;
for i = 1:pars.r
  J = J + ( meta.S(i)/( meta.S(i)^2 + pars.alpha ) ) * ...
    reshape( meta.V(:,i), pars.n, 1 ) * ( meta.U(:,i)' * result.data.Y );
end
solution.algTime = toc(parTic);
solution.J = diag(meta.ColumnNorm.^-1) * J;

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
function U = Tikhonov_Ucurve( meta, result, pars, alpha)
% CRESO

% solution
J = meta.Leadfield' * pinv( eye(pars.m) + alpha * meta.Leadfield * meta.Leadfield' ) * result.data.Y;

% norm
N = vecnorm( J, 2 )^2;

% residual
R = vecnorm( meta.Leadfield*J - result.data.Y, 2 )^2;

U = 1/R + 1/N;

end
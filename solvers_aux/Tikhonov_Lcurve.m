function [N, R] = Tikhonov_Lcurve( meta, result, pars, alpha)
% L-Curve Criterion

% solution
J = meta.LeadfieldOG' * pinv( eye(pars.m) + alpha * meta.LeadfieldOG * meta.LeadfieldOG' ) * result.data.Y;

% norm
N = vecnorm( J, 2 )^2;

% residual
R = vecnorm( meta.Leadfield*J - result.data.Y, 2 )^2;

end
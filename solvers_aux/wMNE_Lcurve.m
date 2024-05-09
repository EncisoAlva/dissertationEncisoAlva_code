function [N, R] = wMNE_Lcurve( meta, result, pars, alpha)
% L-Curve Criterion

% solution
J = meta.Leadfield' * pinv( eye(pars.m) * alpha + meta.Leadfield * meta.Leadfield' ) * result.data.Y;

% norm
N = vecnorm( J, 2 )^2;

% residual
R = vecnorm( meta.Leadfield*J - result.data.Y, 2 )^2;

end
function C = wMNE_CRESO( meta, result, pars, alpha)
% CRESO

% solution
J = meta.Leadfield' * pinv( eye(pars.m) * alpha + meta.Leadfield * meta.Leadfield' ) * result.data.Y;

% norm
N = vecnorm( J, 2 )^2;

% residual
R = vecnorm( meta.Leadfield*J - result.data.Y, 2 )^2;

C = -R + alpha*N;
%C = -R + N;

end
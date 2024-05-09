function G = wMNE_GCV( meta, result, pars, alpha)
% L-Curve Criterion

% kernel
K = meta.Leadfield' * pinv( eye(pars.m) * alpha + meta.Leadfield * meta.Leadfield' );

% solution
J = K * result.data.Y;

% residual
R = vecnorm( meta.Leadfield*J - result.data.Y, 2 )^2;

% trace
T = ( pars.m - trace( meta.Leadfield* K ) )^2;

% GCV
G = pars.m * R / T;

end
function G = wMNE_GCV( meta, result, pars, alpha)
% L-Curve Criterion

% kernel
K = meta.LeadfieldColNorm' * pinv( eye(pars.m) * alpha + meta.LeadfieldColNorm * meta.LeadfieldColNorm' );

% solution
J = K * (result.data.Y) ./(meta.ColumnNorm');

% residual
R = vecnorm( meta.LeadfieldColNorm*J - result.data.Y, 2 )^2;

% trace
T = ( pars.m - trace( meta.LeadfieldColNorm* K ) )^2;

% GCV
G = pars.m * R / T;

end
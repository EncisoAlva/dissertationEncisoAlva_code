function U = sLORETA_Ucurve( meta, info, result, pars, alpha )
% Generalized Cross-Validation

% inversion kernel
Kv = meta.LeadfieldColNorm' * pars.H * pinv( pars.HGGH + alpha*(pars.H) );

% solution
J = Kv * result.data.Y;

% norm
N = vecnorm( J, 2 )^2;

% residual
R = vecnorm( meta.LeadfieldColNorm*J - result.data.Y, 2 )^2;

% U-curve metric
U = 1/R + 1/(alpha* N);

end
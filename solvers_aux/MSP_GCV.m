function G = MSP_GCV( meta, info, result, pars, alpha )
% Generalized Cross-Validation

% inversion kernel
Kv = diag( pars.W.^(-2) ) * meta.LeadfieldColNorm' * pinv( pars.GWG + alpha*eye(pars.m) );

% GCV metric
G = pars.m * norm( (meta.LeadfieldColNorm*Kv-eye(pars.m)) *result.data.Y,"fro")^2 ...
  / ( trace(meta.LeadfieldColNorm*Kv) - pars.m )^2;
end
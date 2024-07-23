function G = sLORETA_GCV( meta, info, result, pars, alpha )
% Generalized Cross-Validation

% inversion kernel
Kv = meta.LeadfieldColNorm' * pars.H * pinv( pars.HGGH + alpha*pars.H );

% GCV metric
G = norm( meta.LeadfieldColNorm*Kv*result.data.Y - ...
  result.data.Y,"fro")^(1/2) ...
  / ( trace(meta.LeadfieldColNorm*Kv) - pars.M )^2;
end
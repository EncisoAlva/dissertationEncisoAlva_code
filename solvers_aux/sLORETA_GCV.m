function G = sLORETA_GCV( meta, info, result, pars, alpha )
% Generalized Cross-Validation

% inversion kernel
Kv = meta.Leadfield' * pars.H * pinv( pars.HGGH + alpha*eye(pars.M) );

% GCV metric
G = pars.M * norm( (meta.Leadfield*Kv - pars.H ) * ...
  result.data.Y,"fro")^2 ...
  / ( trace(meta.Leadfield*Kv) - pars.M )^2;
end
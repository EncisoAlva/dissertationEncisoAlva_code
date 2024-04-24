function G = SingleRegionPrior_GCV( meta, info, result, pars, lambda, Hs )
% Generalized Cross-Validation

% inversion kernel
Kv = Hs * meta.Leadfield' * ...
    pinv( meta.Leadfield*Hs*meta.Leadfield' + lambda*eye(pars.M) ) ;

% GCV metric
G = pars.M * norm( (meta.Leadfield*Kv-eye(pars.M)) *result.data.Y,"fro")^2 ...
  / ( trace(meta.Leadfield*Kv) - pars.M )^2;
end
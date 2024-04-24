function solution = sLORETA( meta, result, buffer )
  % common parameters
  meta.M = size(meta.Leadfield, 1);
  meta.N = size(meta.Leadfield, 2);
  meta.T = size(result.data.time, 2);

  meta.H = eye(meta.N) - ones(meta.N)/meta.N;

  % hyperparameter tuning via Generalized Cross-Validation
  parTic  = tic;
  params  = sLORETA_optimize(meta, result);
  parTime = toc(parTic);
  % solution per se
  algTic   = tic;
  solution = sLORETA_solve(   meta, result, params);
  algTime  = toc(algTic);
  %
  solution.parTime = parTime;
  solution.algTime = algTime;
  %solution.svdTime = svdTime;
  solution.params = params;

  solution.buffer = [];
end

function solution = sLORETA_solve( meta, result, params )
  kernel = sLORETA_kernel( meta, params.alpha );

  solution.J = kernel* result.data.Y_OG;
end

function params =  sLORETA_optimize( meta, result)
  params = [];

  alphas = (median(meta.S)^2) * ( 10.^(-10:10) );
  Gs     = zeros( size(alphas) );

  found = false; iter = 0; big = 1;
  while ~found && (iter<2)
    for q = 1:length(alphas)
      alpha = alphas(q);
      Gs(q)  = sLORETA_GCV( meta, result, alpha );
    end
    [~, idx] = min(Gs);
    if (idx==1) || (idx==21)
      alphas = (alphas(idx)) * ( 10.^(-10:10) );
      Gs     = zeros( size(alphas) );
    else
      alphas = (alphas(idx)) * ( 10.^((-big):(big/10):big) );
      big    = big/10;
      iter   = iter + 1;
    end
  end
  [~, idx] = min(Gs);
  params.alpha = alphas(idx);
end

function G = sLORETA_GCV( meta, result, alpha)
  Kv = sLORETA_kernel( meta, alpha );
  G  = meta.M*norm( (meta.Leadfield*Kv-eye(meta.M)) *result.data.Y,"fro")^2 ...
    / trace(meta.Leadfield*Kv-eye(meta.M))^2;
end

function kernel = sLORETA_kernel( meta, alpha )
%
  kernel = (meta.Leadfield') * pinv( (meta.Leadfield)*(meta.Leadfield') + alpha*eye(meta.M) );
end
function solution = Tikhonov( meta, result, buffer )
  % common parameters
  %tic
  %[U,S,V] = svd(meta.Leadfield);
  %svdTime = toc;
  %meta.U  = U;
  %meta.S  = diag(S);
  %meta.V  = V;
  %
  meta.m = size(meta.Leadfield, 1);
  meta.n = size(meta.Leadfield, 2);
  meta.r = min(meta.m, meta.n);
  meta.t = size(result.data.time, 2);

  % hyperparameter tuning via Generalized Cross-Validation
  parTic  = tic;
  params  = Tikhonov_optimize(meta, result);
  parTime = toc(parTic);
  % solution per se
  algTic   = tic;
  solution = Tikhonov_solve(   meta, result, params);
  algTime  = toc(algTic);
  %
  solution.parTime = parTime;
  solution.algTime = algTime;
  %solution.svdTime = svdTime;
  solution.params = params;

  solution.buffer = [];
end

function solution = Tikhonov_solve( meta, result, params )
  solution = [];
  J = zeros( meta.n, meta.t );
  for i = 1:meta.r
    J = J + ( meta.S(i)/( meta.S(i)^2 + params.alpha ) ) * ...
    reshape( meta.V(:,i), meta.n, 1 ) * ( meta.U(:,i)' * result.data.Y_OG );
  end
  solution.J = J;
end

function params =  Tikhonov_optimize( meta, result)
  params = [];

  alphas = (median(meta.S)^2) * ( 10.^(-10:10) );
  Gs     = zeros( size(alphas) );

  found = false; iter = 0; big = 1;
  while ~found && (iter<2)
    for q = 1:length(alphas)
      alpha = alphas(q);
      Gs(q)  = Tikhonov_GCV( meta, result, alpha );
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

function G = Tikhonov_GCV( meta, result, alpha)
  G = 0;
  for i = 1:meta.r
    G = G + mean( ( (alpha/(meta.S(i)^2+alpha)) *  meta.U(:,i)' * result.data.Y ).^2, 2);
  end
  for i = (meta.r+1):meta.m
    G = G + sum( ( meta.U(:,i)' * result.data.Y ).^2 );
  end
  den = 0;
  for i = 1:meta.r
    den = den + (( meta.S(i)^2 )/( meta.S(i)^2 + alpha ));
  end
  G = G /( (meta.m - den)^2 );
end
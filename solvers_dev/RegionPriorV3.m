function solution = RegionPriorV3( meta, result, buffer )
% this path simplifies to only 1 region

  % common parameters
  meta.M = size(meta.Leadfield, 1);
  meta.N = size(meta.Leadfield, 2);
  meta.T = size(result.data.time, 2);

  %if meta.optimizeParams
  %meta.K = size( meta.atlas,2 ) +1;
  meta.K = 2;

  % Region activations
  meta.SS = false(meta.K,1);
  meta.n  = zeros(meta.K,1);
  if meta.constrained
    meta.R = cell(   (meta.K),1 );
  else
    meta.R = cell( 3*(meta.K),1 );
  end
  if meta.constrained
    IDX = 1:(meta.N);
    meta.R{1} = IDX(vecnorm( meta.Gridloc - result.data.SCtr, 2, 2) <  3*result.data.tau);
    meta.R{2} = IDX(vecnorm( meta.Gridloc - result.data.SCtr, 2, 2) >= 3*result.data.tau);
    %
    meta.n(1) = size(meta.R{1},2);
    meta.n(2) = size(meta.R{2},2);
  else
    IDX = 1:(meta.N);
    for nu = 1:3
      meta.R{  nu} = IDX(vecnorm( meta.Gridloc - result.data.SCtr, 2, 2) <  3*result.data.tau)*3-3+nu;
      meta.R{3+nu} = IDX(vecnorm( meta.Gridloc - result.data.SCtr, 2, 2) >= 3*result.data.tau)*3-3+nu;
      %
      meta.n(  nu) = size(meta.R{  nu},2);
      meta.n(3+nu) = size(meta.R{3+nu},2);
    end
  end
  %unassigned = true(meta.N,1);
  %for k = 1:(meta.K-1)
  %  if meta.constrained
  %    meta.R{k} = meta.atlas(k).Vertices;
  %    unassigned(meta.R{k}) = false;
  %    meta.n(k) = size(meta.R{k},2);
  %  else
  %    for nu = 1:3
  %      meta.R{3*k-3+nu} = meta.atlas(k).Vertices*3-3+nu;
  %      unassigned(meta.R{3*k-3+nu}) = false;
  %    end
  %    meta.n(k) = size(meta.R{3*k-1},2);
  %  end
  %end
  %if any(unassigned)
  %  idx = 1:(meta.N);
  %  if meta.constrained
  %    idx = unique(idx(unassigned));
  %    meta.R{meta.K} = idx;
  %    meta.n(meta.K) = size(meta.R{meta.K},2);
  %  else
  %    idx = unique(ceil(idx(unassigned)/3-1));
  %    for nu = 1:3
  %      meta.R{3*meta.K-3+nu} = idx*3-3+nu;
  %    end
  %    meta.n(meta.K) = size(meta.R{meta.K*3-1},2);
  %  end
  %else
  %  meta.K = meta.K-1;
  %  meta.n(meta.K) = [];
  %end

  %
  tmp = [];
  tmp.buffer.M  = meta.M;
  tmp.buffer.N  = meta.N;
  tmp.buffer.T  = meta.T;
  tmp.buffer.K  = meta.K;
  %tmp.buffer.SS = meta.SS;
  %tmp.buffer.n  = meta.n;
  %tmp.buffer.R  = meta.R;
  %else
  %  meta.M  = buffer.M;
  %  meta.N  = buffer.N;
  %  meta.T  = buffer.T;
  %  meta.K  = buffer.K;
  %  meta.SS = buffer.SS;
  %  meta.n  = buffer.n;
  %  meta.R  = buffer.R;
  %end

  % using knowledge of the true active region
  %meta.SS(result.meta.idxRegion) = true;
  meta.SS(1) = true;
  
  %
  if meta.constrained
    meta.Ls = zeros(meta.N, meta.K);
    for k = 1:meta.K
      if meta.SS(k)
        meta.Ls(meta.R{k},k) = 1;
      end
    end
  else
    meta.Ls = zeros(meta.N, 3*meta.K);
    for k = 1:meta.K
      if meta.SS(k)
        for nu = 1:3
          meta.Ls(3*k-3+nu,3*k-3+nu) = 1;
        end
      end
    end
  end

  meta.constrained = false;

  % hyperparameter tuning via Generalized Cross-Validation
  parTic  = tic;
  if meta.optimizeParams
    %params = [];
    %params.gamma0 = 1;
    %params.gamma1 = 1;
    params  = RegionPrior_optimize(meta, result);
  else
    params = buffer.params;
  end
  parTime = toc(parTic);
  % solution per se
  algTic   = tic;
  solution = RegionPrior_solve(   meta, result, params);
  algTime  = toc(algTic);
  %
  solution.parTime = parTime;
  solution.algTime = algTime;
  %solution.svdTime = svdTime;
  solution.params = params;
  if meta.optimizeParams
    solution.buffer        = tmp.buffer;
    solution.buffer.params = params;
  end
end

function solution = RegionPrior_solve( meta, result, params )
  kernel = RegionPrior_kernel( meta, params );

  solution.J = kernel* result.data.Y_OG;
end

function params =  RegionPrior_optimize( meta, result)
  params = [];
  par = [];

  gammas = [ (median(meta.S)^2) * ( 10.^(-10:10) ) ; (median(meta.S)^2) * ( 10.^(-10:10) )];
  Gs     = zeros( size(gammas,2), size(gammas,2) );

  found = false; iter = 0; 
  big = [1,1];
  while ~found && (iter<3)
    for q0 = 1:size(gammas,2)
    for q1 = 1:size(gammas,2)
      par.gamma0 = gammas(1,q0);
      par.gamma1 = gammas(2,q1);
      Gs(q0,q1)  = RegionPrior_GCV( meta, result, par );
    end
    end
    [~, IDX] = min(Gs,[],"all");
    [idx1, idx2 ] = ind2sub( size(Gs) ,IDX);
    idx = [idx1, idx2 ];
    for qq = 1:2
    if (1 == idx(qq)) || (idx(qq) == 21)
      gammas(qq,:) = gammas(qq,idx(qq)) * ( 10.^(-big(qq):(big(qq)/10):big(qq)) );
    else
      gammas(qq,:) = gammas(qq,idx(qq)) * ( 10.^(-big(qq):(big(qq)/10):big(qq)) );
      big(qq) = big(qq)/10;
      iter   = iter + 1/2;
    end
    end
    %gammas0 = (gammas0(idx)) * ( 10.^((-big):(big/10):big) );
    %found=true;
  end
  [~, IDX] = min(Gs,[],"all");
  [idx1, idx2 ] = ind2sub( size(Gs) ,IDX);

  params.gamma0 = gammas(1,idx1);
  params.gamma1 = gammas(2,idx2);
end

function G = RegionPrior_GCV( meta, result, par)
  Kv = RegionPrior_kernel( meta, par );
  G  = meta.M*norm( (meta.Leadfield*Kv-eye(meta.M)) *result.data.Y,"fro")^2 ...
    / trace(meta.Leadfield*Kv-eye(meta.M))^2;
end

function kernel = RegionPrior_kernel( meta, par )
  if meta.constrained
    tmpK =   meta.K;
  else
    tmpK = 3*meta.K;
  end
  LLinv = zeros(tmpK, tmpK);
  for k = 1:tmpK
    LLinv(k,k) = 1/meta.n(k);
  end
  key    = ( par.gamma1* meta.Leadfield * meta.Ls*LLinv*meta.Ls' + par.gamma0*meta.Leadfield ) * meta.Leadfield' ...
    + eye(meta.M);
  keyInv = linsolve( key, eye(meta.M) );
  kernel = (par.gamma0 * meta.Leadfield' + par.gamma1* meta.Ls*LLinv*meta.Ls' * meta.Leadfield') ...
    * keyInv;
end
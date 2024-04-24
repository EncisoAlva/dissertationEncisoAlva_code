function pars = zSISSY_tune( meta, info, result )
% TODO add description (optional)
%
% J^ = armin_J || G*J - Y ||^2_F + lambda*( || V*J ||_1 + alpha*|| J ||_1 )

% start timer
parTic = tic;
pars   = [];

% general values
pars.m  = size(meta.Leadfield, 1);
pars.n  = size(meta.Leadfield, 2);
pars.t  = size(result.data.time, 2);

% total variation operator
switch info.SourceType
  case 'surface'
    pars.Edge = edges_mine( meta.Cortex );
  case 'volume'
    pars.Edge = edges_mine( meta. DT );
end
pars.nEdges = size(pars.Edge,1);
%meta.V = sparse( pars.nEdges, pars.n );
%for ii = 1:size(pars.Edge,1)
%  V( ii, pars.Edge(ii,1) ) =  1;
%  V( ii, pars.Edge(ii,2) ) = -1;
%end
pars.V = sparse( [(1:pars.nEdges)' ; (1:pars.nEdges)'], ...
  [ pars.Edge(:,1); pars.Edge(:,2)], ...
    [ones(pars.nEdges,1); -ones(pars.nEdges,1)], ...
    size(pars.Edge,1), pars.n );

% parameters recommended on the original paper
pars.alpha = 0.07;
pars.rho   = 1;

% performing an inversion one single time
%pars.invP = pinv( meta.Leadfield' * meta.Leadfield + ...
%    pars.rho*( pars.V'*pars.V + sparse(eye( pars.n )) )  );

% using Woodbury inversion lemma
pars.MM    = sparse( pars.rho*( pars.V'*pars.V + sparse(eye( pars.n )) ) );
pars.invMM = sparse( inv(pars.MM) );
pars.invP  = pars.invMM - pars.invMM * meta.Leadfield' * ...
  pinv( eye(pars.m) + meta.Leadfield * pars.invMM * meta.Leadfield' ) * ...
  meta.Leadfield * pars.invMM;

% hyperparameter tuning by L-curve criterion
best_lambda  = median(meta.S)^2;

% hot start: the MNE solution
admm = [];
%admm.J = meta.Leadfield' * ...
%  pinv( meta.Leadfield * meta.Leadfield' + best_lambda* eye(pars.m) ) * result.data.Y;
%admm.J(abs(admm.J) < 0.1*max(abs(admm.J))) = 0;
%admm.J = admm.J / max(abs(admm.J));

% cold start: zeros
admm.J = zeros(pars.n, 1);
admm.A = pars.V * admm.J;
admm.B = admm.J;
admm.C = zeros(size(admm.A));
admm.D = zeros(size(admm.B));
admm.iter = 0;
admm.res  = 0;

if info.debugFigs
  figure()
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
    40, abs(admm.J)*120,'filled')
  colormap("parula")
   scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
    200, 'red','filled')
  title("Magnitude of true sources")
  b = colorbar;
  b.Label.String = 'Unitless; range=[0,120]';
end

% search
for iter = 1:1
  % try many values for lambda
  lambdas = best_lambda * (2.^( linspace(-30,30,60) ));
  Gs      = zeros( size(lambdas) );
  for q = 1:length(lambdas)
    lambda  = lambdas(q);
    % if the previous iteration was bad, re-initialize
    if admm.res > 10^10
      admm.J = meta.Leadfield' * ...
        pinv( meta.Leadfield * meta.Leadfield' + best_lambda* eye(pars.m) ) * result.data.Y;
      admm.A = pars.V * admm.J;
      admm.B = admm.J;
      admm.C = zeros(size(admm.A));
      admm.D = zeros(size(admm.B));
      admm.iter = 0;
    end
    admm  = SISSY_ADMM( meta, result, pars, pars.alpha, lambda, admm );
    %
    % J^ = armin_J || G*J - Y ||^2_F + lambda*( || V*J ||_1 + alpha*|| J ||_1 )
    Gs(q) = norm([ norm(meta.Leadfield*admm.J-result.data.Y,2), ...
      norm(pars.V*admm.J,1) + pars.alpha*norm(admm.J,1) ], 1);
  end
  [~, idx]   = min(Gs);
  best_lambda = lambdas(idx);
end
pars.lambda   = max(best_lambda, 0.001);
% for a future hot-start
pars.wMNE_ker = meta.Leadfield' * ...
  pinv( meta.Leadfield * meta.Leadfield' + best_lambda* eye(pars.m) );

% report total number of iterations in tuning
pars.tuning_iterations = admm.iter;

% stop timer
pars.parTime = toc(parTic);

end
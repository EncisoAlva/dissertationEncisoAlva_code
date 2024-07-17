function solution = zSISSY( meta, info, result, pars )
% TODO add description (optional)
%
% J^ = armin_J || G*J - Y ||^2_F + lambda*( || V*J ||_1 + alpha*|| J ||_1 )
%
solution = [];
parTic = tic;

% hot start
admm = [];
admm.J = pars.wMNE_ker * result.data.Y;
admm.J(abs(admm.J) < 0.1*max(abs(admm.J))) = 0;
admm.A = pars.V * admm.J;
admm.B = admm.J;
admm.C = zeros(size(admm.A));
admm.D = zeros(size(admm.B));
admm.iter = 0;

%solution using ADMM
admm = SISSY_ADMM( meta, result, pars, pars.alpha, pars.lambda, admm );
J    = admm.J;
solution.J = diag(meta.ColumnNorm.^-1) * J;
solution.iter = admm.iter;

% norm of J
switch meta.Type
  case 'surface'
    solution.normJ = abs( solution.J );
  case 'volume'
    solution.normJ = dip_norm( solution.J );
end

% stop timer
solution.algTime = toc(parTic);
end
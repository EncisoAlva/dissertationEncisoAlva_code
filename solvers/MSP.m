function solution = MSP( meta, info, result, pars )
% TODO add description (optional)
%
solution = [];

% solution per se
parTic = tic;

% project measurements over the first s eigenvectors of G
Us = meta.U(:, (1:pars.s));
Ys = Us * Us' * result.data.Y / norm( result.data.Y );

% project G into the projection of Y
Ps = Ys * pinv( Ys'*Ys ) * Ys';

% Activation Probability Map (APM) is defined as follows
D = zeros( meta.nGridDips, 1 );
for ii = 1:meta.nGridDips
  D(ii) = norm( Ps * meta.Leadfield(:,ii), 2 )^2;
end
APM = D;

% another thing to do with APM is to create weights
W   = 1 - APM *.95; % robustness(?)
GWG = meta.Leadfield * diag( W.^(-2) ) * meta.Leadfield';

% solution per se
kernel = diag( W.^(-2) ) * meta.Leadfield' * pinv( GWG + pars.alpha*eye(pars.m) );
J = kernel* result.data.Y;
solution.algTime = toc(parTic);

solution.J = diag(meta.ColumnNorm.^-1) * J;

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
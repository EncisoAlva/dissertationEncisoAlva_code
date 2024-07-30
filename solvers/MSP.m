function solution = MSP( meta, info, result, pars )
% TODO add description (optional)
%
solution = [];

% solution per se
parTic = tic;

% project measurements over the first s eigenvectors of G
Us = pars.Unorm(:, (1:pars.s));
Ys = Us * Us' * result.data.Y / norm( result.data.Y );

% project G into the projection of Y
Ps = Ys * pinv( Ys'*Ys ) * Ys';

% Activation Probability Map (APM) is defined as follows
D = zeros( meta.nGridDips, 1 );
switch info.SourceType
  case 'surface'
    for ii = 1:meta.nGridDips
      D(ii) = norm( Ps * meta.LeadfieldColNorm(:,ii), 2 )^2;
    end
  case 'volume'
    for ii = 1:meta.nGridDips
      D(ii) = 0;
      for tt = 1:3
          D(ii) = D(ii) + norm( Ps * meta.LeadfieldColNorm(:,3*(ii-1)+tt), 2 )^2/3;
      end
    end
end
APM = D;

% another thing to do with APM is to create weights
W   = 1 - APM *.99; % robustness(?)
switch info.SourceType
  case 'surface'
    % nothing else
  case 'volume'
    W = kron( W, [1,1,1]' );
end
GWG = meta.LeadfieldColNorm * diag( W.^(-2) ) * meta.LeadfieldColNorm';

% solution per se
kernel = diag( W.^(-2) ) * meta.LeadfieldColNorm' * pinv( GWG + pars.alpha*eye(pars.m) );
J = kernel* result.data.Y;
solution.algTime = toc(parTic);

%solution.J = diag(meta.ColumnNorm.^-1) * J;
solution.J = J;

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
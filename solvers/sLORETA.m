function solution = sLORETA( meta, info, result, pars )
% TODO add description (optional)
%
solution = [];

% solution per se
parTic = tic;
J = pars.kernel* result.data.Y;
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
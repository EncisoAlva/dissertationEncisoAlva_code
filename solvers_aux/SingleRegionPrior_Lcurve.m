function metric = SingleRegionPrior_Lcurve( meta, info, result, pars, ...
  lambda, Hs, g0, g1, As )
% L-curve criterion, based on the optimization formulation
%      J^ = B^ + L * R^
%   B^,R^ = argmin || G*(B + L*R) - Y ||^2 + LAM*( ||B||^2/g0 + ||L*R||^2/g1 ) 

% solution depending on lambda
J = As * meta.Leadfield' * ...
    pinv( meta.Leadfield*Hs*meta.Leadfield' + lambda*eye(pars.M) ) *result.data.Y;

% norm of the residual, || G*J - Y ||^2
resNorm = norm( meta.Leadfield * J - result.data.Y, 'fro' );

% norm of J = B + L*R, with 
%   B = H*J      background activity
%   R = L*L'/Nk  regional average
B  = As * J;
LR = J - B;

varNorm = norm( B, 'fro' )^2/g0 + norm( LR, 'fro' )^2/g1;

% L-curve metric: distance to origin
metric = norm( [ resNorm, varNorm ] );
end
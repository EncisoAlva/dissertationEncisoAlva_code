function G = Tikhonov_GCV( meta, result, pars, alpha)
% Generalized Cross-Validation given the SVD decomposition

% residual of Y from leave-one-out
G = 0;
for i = 1:pars.r
  G = G + ( (alpha/(meta.S(i)^2+alpha)) *  meta.U(:,i)' * result.data.Y ).^2;
end
for i = (pars.r+1):pars.m
  G = G + sum( ( meta.U(:,i)' * result.data.Y ).^2 );
end
% trace
tra = 0;
for i = 1:pars.r
  tra = tra + (( meta.S(i)^2 )/( meta.S(i)^2 + alpha ));
end
G = G /( (pars.m - tra)^2 );

end
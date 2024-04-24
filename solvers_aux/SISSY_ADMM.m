function admm_new = SISSY_ADMM( meta, result, pars, alpha, lambda, admm_old)
% solving the SISSY problem with ADMM

% hot start is provided
J_old = admm_old.J;
A_old = admm_old.A;
B_old = admm_old.B;
C_old = admm_old.C;
D_old = admm_old.D;

% loop
residual = 1; counter = 0;
while (residual > 10^-7)&&(residual < 10^15)&&( counter < 2000 ) % arbitrary tolerance
  % update step
  J_new = pars.invP*( meta.Leadfield'*result.data.Y + ...
      pars.rho*pars.V'*(A_old-C_old) + pars.rho*(B_old-D_old) );
  A_new = shrink( lambda/pars.rho,       pars.V*J_new + C_old );
  B_new = shrink( lambda*alpha/pars.rho,        J_new + D_old );
  C_new = C_old + pars.V*J_new - A_new;
  D_new = D_old +        J_new - B_new;
  % report error and move the variables
  %err = norm( J_new-J_old, 1 ) / norm( J_old, 1 );
  residual = norm( meta.Leadfield*(J_new-J_old), 2 );
  maxJ = max(abs(J_new));
  fprintf("Iter %3d. Err = %3.3d. (0.00): %5d dips. (0.05) : %5d dips.  (0.50) : %5d dips.\n", ...
    counter, residual, sum(J_new~=0), sum(abs(J_new)>0.05*maxJ), sum(abs(J_new)>0.5*maxJ))
  %histogram(J_new)
  %
  J_old = J_new;
  A_old = A_new;
  B_old = B_new;
  C_old = C_new;
  D_old = D_new;
  counter = counter + 1;
  %
  if false
    figure()
     trisurf(meta.Cortex.Faces, ...
      meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
    hold on
    scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
      40, abs(J_new)*120,'filled')
    colormap("parula")
    scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
      200, 'red','filled')
    title("Magnitude of true sources")
    b = colorbar;
    b.Label.String = 'Unitless; range=[0,120]';
  end
end

% saving the results
admm_new = [];
admm_new.J = J_new;
admm_new.A = A_new;
admm_new.B = B_new;
admm_new.C = C_new;
admm_new.D = D_new;
admm_new.iter = admm_old.iter + counter;
admm_new.res  = residual;
end
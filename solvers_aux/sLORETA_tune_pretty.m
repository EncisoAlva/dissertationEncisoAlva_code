function pars = sLORETA_tune_pretty( meta, info, result )
% TODO add description (optional)
%

% start timer
parTic = tic;
pars   = [];

% general values
pars.M    = size(meta.Leadfield, 1);
pars.N    = size(meta.Leadfield, 2);
pars.T    = size(result.data.time, 2);
pars.H    = eye(pars.M) - (ones(pars.M)/pars.M); % averaging operator
pars.HGGH = pars.H * meta.Leadfield * meta.Leadfield' * pars.H;

%%
% hyperparameter tuning via L-curve
alphas = 10.^(1:0.1:7);
Rs     = zeros( size(alphas) );
Ns     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  [Ns(q), Rs(q)] = sLORETA_Lcurve( meta, info, result, pars, alpha );
end

figure()
tiledlayout(2,2,'Padding','tight');
nexttile
loglog(Rs,Ns)
grid on
xlabel('$\log \Vert \mathbf{G}\,\mathbf{J}_\lambda - \mathbf{Y} \Vert_2^2 $','Interpreter','latex')
ylabel('$\log \Vert \mathbf{J}_\lambda \Vert_2^2$','Interpreter','latex')
hold on
for aa = 10.^(1:1:7)
  [~,qidx] = min( abs(alphas-aa) );
  text(Rs(qidx), Ns(qidx),  ...
    ['$\lambda = 10^{',num2str(log10(aa)),'}$'], 'Interpreter','latex',...
    'BackgroundColor', 'white')
  scatter(Rs(qidx), Ns(qidx), 30,'blue','filled')
end
title('L-curve')

% hyperparameter tuning via Generalized Cross-Validation
alphas = 10.^(1:0.1:7);
Gs     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Gs(q) = sLORETA_GCV( meta, info, result, pars, alpha );
end

[~, idx]   = min(Gs);
best_alpha = alphas(idx);

nexttile
loglog(alphas,Gs)
grid on
xlabel('$\log \lambda$','Interpreter','latex')
%ylabel('$\log {GCV}(\lambda)$','Interpreter','latex')
ylabel('$\log( \Vert \mathbf{G}\, \mathbf{J}_\lambda - \mathbf{Y} \Vert_F\, /\, tr( \mathbf{G}\, \mathbf{K}_\lambda - \mathbf{I}_M ) )^2$','Interpreter','latex')
hold on
scatter(best_alpha, Gs(idx),50,'red','filled','pentagram')
title('GCV')
xlim([min(alphas) max(alphas)])

% hyperparameter tuning via CRESO
%alphas = 10.^(-4:0.05:5);
Cs     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Cs(q) = sLORETA_CRESO( meta, info, result, pars, alpha );
end

dCs = zeros( size(alphas) );
dCs(1)   = ( Cs(2)   - Cs(1)     )/(alphas(2)   - alphas(1));
dCs(end) = ( Cs(end) - Cs(end-1) )/(alphas(end) - alphas(end-1));
for q = 2:(length(alphas)-1)
  dCs(q) = ( Cs(q+1) - Cs(q-1) )/(alphas(q+1) - alphas(q-1));
end

idx = find( dCs>0, 1, 'last' );
%[~, idx]   = max(Cs);
best_alpha = alphas(idx);

nexttile
semilogx(alphas(2:(end-1)),dCs(2:(end-1)))
hold on
grid on
xlabel('$\log \lambda$','Interpreter','latex')
ylabel('$\frac{d}{d\lambda} ( -\Vert \mathbf{G}\, \mathbf{J}_\lambda - \mathbf{Y} \Vert_2^2 - \lambda \Vert \mathbf{J}_\lambda \Vert_2^2)$','Interpreter','latex')
scatter(best_alpha, dCs(idx),50,'red','filled','pentagram')
title('CRESO')
xlim([min(alphas) max(alphas)])

% hyperparameter tuning via U-curve
%alphas = 10.^(-10:0.1:10);
Us     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Us(q) = sLORETA_Ucurve( meta, info, result, pars, alpha );
end

[~, idx]   = min(Us);
best_alpha = alphas(idx);

nexttile
loglog(alphas,Us)
grid on
xlabel('$\log \lambda$','Interpreter','latex')
ylabel('$\log ( \Vert \mathbf{G}\, \mathbf{J}_\lambda - \mathbf{Y} \Vert_2^{-2} + \Vert \mathbf{J}_\lambda \Vert_2^{-2} )$','Interpreter','latex')
hold on
scatter(best_alpha, Us(idx),50,'red','filled','pentagram')
title('U-curve')
xlim([min(alphas) max(alphas)])

% save graph
set(gcf,'color','w');
fig = gcf;
fig.Units = 'inches';
fig.OuterPosition = [0 0 6 6];

exportgraphics(fig,'ParTuning_sLORETA.pdf','Resolution',300)

%%
% hyperparameter tuning via Generalized Cross-Validation
% starting at median eigenvalue
best_alpha = median(meta.S)^2;
scale  = 10;
for iter = 1:6
  % try many values for alpha, compute GCV value for each, get the min
  alphas = best_alpha * (2.^( (-scale):(scale/10):scale ));
  Gs     = zeros( size(alphas) );
  for q = 1:length(alphas)
    alpha = alphas(q);
    Gs(q) = sLORETA_GCV( meta, info, result, pars, alpha );
  end
  [~, idx]   = min(Gs);
  best_alpha = alphas(idx);
  %
  % if not on the border, reduce scale; else, increase it
  if (1<idx) && (idx<length(Gs))
    scale = scale/10;
  else
    scale = scale*10;
  end
end
pars.alpha  = max(best_alpha, 0.001);
pars.kernel = meta.Leadfield' * pars.H * pinv( pars.HGGH + pars.alpha*eye(pars.M) );

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for sLORETA solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
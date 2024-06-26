function pars = wMNE_tune_pretty( meta, info, result )
% Tuning via Generalized Cross-Validation for wMNE, and additional
% initialization routines.
%
% Weighten Minimum-Norm Estimator (wMNE), follows the basic Tikonov
% regularized estimation
%   J^ = argmin_J || G*J-Y ||^2_F + alpha || W*J ||^2_F
% with W the weight induced by column-normalization of G and ||*||_F is the
% Frobenius norm.
%

% start timer
parTic = tic;
pars   = [];

% general values
pars.m = size(meta.LeadfieldColNorm, 1);
pars.n = size(meta.LeadfieldColNorm, 2);
pars.r = min(pars.m, pars.n);
pars.t = size(result.data.time, 2);

%% hyperparameter tuning via L-curve criterion

alphas = 10.^(-10:0.1:3);
Rs     = zeros( size(alphas) );
Ns     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  [Nq, Rq] = wMNE_Lcurve( meta, result, pars, alpha );
  Ns(q) = Nq;
  Rs(q) = Rq;
  %[Ns(q), Rs(q)] = wMNE_Lcurve( meta, result, pars, alpha );
end

figure()
t = tiledlayout(2,2,'Padding','tight');
nexttile
loglog(Rs,Ns)
grid on
xlabel('$\log \Vert \mathbf{G}\,\mathbf{J}_\lambda - \mathbf{Y} \Vert_2^2 $','Interpreter','latex')
ylabel('$\log \Vert \mathbf{J}_\lambda \Vert_2^2$','Interpreter','latex')
hold on
for aa = 10.^(0:1:3)
  [~,qidx] = min( abs(alphas-aa) );
  text(Rs(qidx), Ns(qidx),  ...
    ['$\lambda = 10^{',num2str(log10(aa)),'}$'], 'Interpreter','latex',...
    'BackgroundColor', 'white')
  scatter(Rs(qidx), Ns(qidx), 30,'blue','filled')
end
title('L-curve')

%% hyperparameter tuning via Generalized Cross-Validation
% starting at median eigenvalue

alphas = 10.^(-5:0.1:5);
Gs     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Gs(q) = wMNE_GCV( meta, result, pars, alpha );
end

[~, idx]   = min(Gs);
best_alpha = alphas(idx);
GCV_alpha  = best_alpha;

nexttile
loglog(alphas,Gs)
grid on
xlabel('$\log \lambda$','Interpreter','latex')
%ylabel('$\log {GCV}(\lambda)$','Interpreter','latex')
ylabel('$\log( \Vert \mathbf{G}\, \mathbf{J}_\lambda - \mathbf{Y} \Vert_F\, /\, tr( \mathbf{G}\, \mathbf{K}_\lambda - \mathbf{I}_M ) )^2$','Interpreter','latex')
hold on
scatter(best_alpha, Gs(idx),50,'red','filled','pentagram')
title('GCV')

%% hyperparameter tuning via CRESO

alphas = 10.^(-4:0.05:5);
Cs     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Cs(q) = wMNE_CRESO( meta, result, pars, alpha );
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


%% hyperparameter tuning via U-curve

alphas = 10.^(-10:0.1:10);
Us     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Us(q) = wMNE_Ucurve( meta, result, pars, alpha );
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


set(gcf,'color','w');
t.Units = 'inches';
t.OuterPosition = [0 0 6 6];

exportgraphics(t,'ParTuning_wMNE.pdf','Resolution',300)

pars.alpha = max(GCV_alpha, 0.001);

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for wMNE solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
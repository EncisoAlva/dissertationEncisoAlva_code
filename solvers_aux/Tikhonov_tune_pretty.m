function pars = Tikhonov_tune_pretty( meta, info, result )
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

% hyperparameter tuning via L-curve criterion

alphas = 10.^(-10:0.1:10);
Rs     = zeros( size(alphas) );
Ns     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  [Nq, Rq] = Tikhonov_Lcurve( meta, result, pars, alpha );
  Ns(q) = Nq;
  Rs(q) = Rq;
  %[Ns(q), Rs(q)] = Tikhonov_Lcurve( meta, result, pars, alpha );
end

figure()
t = tiledlayout(1,2,'Padding','tight');
nexttile
loglog(Rs,Ns)
grid on
xlabel('$\log \lambda$','Interpreter','latex')
ylabel('$\log {GCV}(\lambda)$','Interpreter','latex')
hold on
scatter(best_alpha, Gs(idx),30,'red','filled')

% hyperparameter tuning via Generalized Cross-Validation
% starting at median eigenvalue

alphas = 10.^(-5:0.1:5);
Gs     = zeros( size(alphas) );
for q = 1:length(alphas)
  alpha = alphas(q);
  Gs(q) = Tikhonov_GCV( meta, result, pars, alpha );
end

[~, idx]   = min(Gs);
best_alpha = alphas(idx);

loglog(alphas,Gs)
grid on
xlabel('$\log \lambda$','Interpreter','latex')
ylabel('$\log {GCV}(\lambda)$','Interpreter','latex')
hold on
scatter(best_alpha, Gs(idx),30,'red','filled')
set(gcf,'color','w');
t.Units = 'inches';
t.OuterPosition = [0 0 4 6];

exportgraphics(t,'GCV_wMNE.pdf','Resolution',300)

pars.alpha = max(best_alpha, 0.001);

% stop timer
pars.parTime = toc(parTic);

% print the results nicely
fprintf("Optimization via GCV for wMNE solver.\n Optimal lambda: ")
disp(pars.alpha)
fprintf("\n")

end
%% Test Case: Data generation
n = 500;
m = 100;
d = 20;
b = m/d;
X = zeros(n,m);
W = zeros(m,d);
repprop = 0.25;
S_noise = 1;
X_noise = 0.1;
nsim = 2;
nfolds = 2;
alphas = [0:0.01:0.1,0.2,0.5,1];

X_bysim = cell(nsim,1);
W_bysim = cell(nsim,1);
W_init_bysim  = cell(nsim,1);
Wl_bysim = cell(nsim,1);
Wl_cv_bysim = cell(nsim,1);
best_alpha_bysim = zeros(nsim,1);
criterion_cv_onfulldata_bysim = zeros(nsim,1);
criterion_old_onfulldata_bysim = zeros(nsim,1);

for l = 1:nsim
    % Simulate the dataset
    for i = 1:d
    ind = (1+(i-1)*b):(i*b);
    if rand > 0.25 || i==1
        sig = randn(n,1);
    end
    X(:,ind) = sig*ones(1,b) + 0.1*randn(n,b);
    W(ind,i) = rand(b,1);
    end

    % change W_init here
    %   W_init = rand(m,d);
    %   W_init = Wo+randn(m,d);
        W_init = double(W~=0);
    %   W_init = double(W>0.25);

    L = X*W;
    S = corr(L);
    S_input = S + S_noise .* randn(d,d);
    
    % Perform Cross-validation
    mean_criterion_byalpha = zeros(length(alphas),1);
    % Cross-validation Partition
    CVO = cvpartition(n,'k',nfolds);
    criterion_byfold = zeros(CVO.NumTestSets,1);

    for k = 1:length(alphas)
        for i = 1:CVO.NumTestSets
            trIdx = CVO.training(i);
            teIdx = CVO.test(i);
            [~,criterion] = estW_smooth(X(trIdx,:),X(teIdx,:),S_input,W_init,alphas(k));
            criterion_byfold(i) = criterion;  
        end
        mean_criterion_byalpha(k) = mean(criterion_byfold);
    
    end
[max_criterion, Idx_max_criteirion] = max(mean_criterion_byalpha);
best_alpha = alphas(Idx_max_criteirion);

fprintf('Best alpha is found to be %0.2f\n', best_alpha);  

% Runing on the full data with best tuning parameter found by CV
[Wl_cv,criterion_cv_onfulldata] = estW_smooth(X,X,S_input,W_init,best_alpha);
[Wl,criterion_old_onfulldata] = estW(X,S_input,W_init);

fprintf('Completed for Simulation %i\n', l); 
X_bysim{l} = X;
W_bysim{l} = W;
W_init_bysim{l}  = W_init;
Wl_bysim{l} = Wl;
Wl_cv_bysim{l} = Wl_cv;
best_alpha_bysim(l) = best_alpha;
criterion_cv_onfulldata_bysim(l) = criterion_old_onfulldata;
criterion_old_onfulldata_bysim(l) = criterion_cv_onfulldata;
disp('Results saved.');
end
  
fname = sprintf('results/cv/simres_all_n%d_m%d_d%d_repprop%.2f_Snoise%.2f_Xnoise%0.2f_nsim%d.mat', n, m, d, repprop, S_noise, X_noise,nsim);
save(fname,'X_bysim','W_bysim','W_init_bysim','Wl_bysim','Wl_cv_bysim','best_alpha_bysim','criterion_cv_onfulldata_bysim','criterion_old_onfulldata_bysim','alphas','nfolds','-v7.3');
disp('Success, All Results saved!');

%% Plotting 
Ll_cv = X * Wl_cv;
Sl_cv = corr(Ll_cv);

Ll = X * Wl;
Sl = corr(Ll);

L_init = X * W_init;
S_init = corr(L_init);

figure;
subplot(2,2,1);
imagesc(S.*~eye(d));
title("Gene-wise corr (L)");
subplot(2,2,2);
imagesc(S_init.*~eye(d));
title("Gene wise corr (L_{init})");
subplot(2,2,3);
imagesc(Sl.*~eye(d));
title("Gene wise corr (L_{learned})");
subplot(2,2,4);
imagesc(Sl_cv.*~eye(d));
title("Gene wise corr (L_{cv})");

figure; plot(zscore(W(W~=0)));
hold on; plot(zscore(Wl(W~=0)));
hold on; plot(zscore(Wl_cv(W~=0)));

figure;
subplot(2,2,1);
imagesc(L);
title("L");
subplot(2,2,2);
imagesc(L_init);
title("L_{init}");
subplot(2,2,3);
imagesc(Ll);
title("L_{learned}");
subplot(2,2,4);
imagesc(Ll_cv);
title("L_{cvlearned}");
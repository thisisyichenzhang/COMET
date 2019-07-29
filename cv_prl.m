%% Run simulation, generate data, Cross-validation, perform regularization method and original method
function cv_prl(n,m,d,repprop,L_noise,X_noise,W_noise,nsim,nfolds)
    alphas = [0:0.01:0.1,0.2,0.5,0.8,1];
    X_bysim = cell(nsim,1);
    X_input_bysim = cell(nsim,1);
    W_bysim = cell(nsim,1);
    W_init_bysim  = cell(nsim,1);
    Wl_bysim = cell(nsim,1);
    Wl_cv_bysim = cell(nsim,1);
    best_alpha_bysim = zeros(nsim,1);
    criterion_cv_onfulldata_bysim = zeros(nsim,1);
    criterion_old_onfulldata_bysim = zeros(nsim,1);
    b = m/d;
    for l = 1:nsim
        %Set the seed to be iterations number at the beginning of each iteration.
        rng(l);
        X = zeros(n,m);
        W = zeros(m,d);
        % Simulate the dataset
        for i = 1:d
        ind = (1+(i-1)*b):(i*b);
        if rand > repprop || i==1
            sig = randn(n,1);
        end
        X(:,ind) = sig*ones(1,b) + 0.1*randn(n,b);
        W(ind,i) = rand(b,1);
        end

        % change W_init here
        %   W_init = rand(m,d);
        %   W_init = W+randn(m,d);
        % !!!WORONG!!! !!LACK BRACKET!! : W_init = double(W + (rand(m,d)>(1-W_noise))~=0);
        W_init = double((W + (rand(m,d)>(1-W_noise)))~=0);
        %   W_init = double(W~=0);
        %   W_init = double(W>0.25);

        L = X*W;
        %S = corr(L);
        L_input = L + L_noise .* randn(n,d);
        S_input = corr(L_input);
        X_input = X + X_noise .* rand(n,m);

        % Perform Cross-validation
        mean_criterion_byalpha = zeros(length(alphas),1);
        % Cross-validation Partition
        CVO = cvpartition(n,'k',nfolds);
        criterion_byfold = zeros(CVO.NumTestSets,1);

        for k = 1:length(alphas)
            for i = 1:CVO.NumTestSets
                trIdx = CVO.training(i);
                teIdx = CVO.test(i);
                [~,criterion] = estW_smooth(X_input(trIdx,:),X_input(teIdx,:),S_input,W_init,alphas(k));
                criterion_byfold(i) = criterion;  
            end
            mean_criterion_byalpha(k) = mean(criterion_byfold);

        end
    [~, Idx_max_criteirion] = max(mean_criterion_byalpha);
    best_alpha = alphas(Idx_max_criteirion);

    fprintf('Best alpha is found to be %0.2f\n', best_alpha);  

    % Runing on the full data with best tuning parameter found by CV
    [Wl_cv,criterion_cv_onfulldata] = estW_smooth_prl(X_input,X_input,S_input,W_init,best_alpha);
    %[Wl,criterion_old_onfulldata] = estW(X_input,S_input,W_init);

    fprintf('Completed for Simulation %i\n', l); 
    X_bysim{l} = X;
    X_input_bysim{l} = X_input;
    W_bysim{l} = W;
    W_init_bysim{l}  = W_init;
    %Wl_bysim{l} = Wl;
    Wl_cv_bysim{l} = Wl_cv;
    best_alpha_bysim(l) = best_alpha;
    criterion_cv_onfulldata_bysim(l) = criterion_cv_onfulldata;
    %criterion_old_onfulldata_bysim(l) = criterion_old_onfulldata;
    disp('Results saved.');
    end

    fname = sprintf('results/cv/prltst_n%d_m%d_d%d_repprop%.2f_Lnoise%.2f_Xnoise%0.2f_Wnoise%0.2f_nsim%d.mat', n, m, d, repprop, L_noise, X_noise, W_noise, nsim);
    save(fname,'X_bysim','X_input_bysim','W_bysim','W_init_bysim','Wl_cv_bysim','best_alpha_bysim','criterion_cv_onfulldata_bysim','alphas','nfolds','-v7.3');
    disp('Success, All Results saved!');

end

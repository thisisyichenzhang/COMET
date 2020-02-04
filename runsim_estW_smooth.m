%% Run simulation, generate data, Cross-validation, perform regularization method and original method
function runsim_estW_smooth(n,m,d,repprop,L_noise,X_noise,W_noise,nsim,nfolds)
    alphas = [0:0.01:0.1,0.2,0.5,0.8,1];
    %X_bysim = cell(nsim,1);
    %X_input_bysim = cell(nsim,1);
    %W_bysim = cell(nsim,1);
    %W_init_bysim  = cell(nsim,1);
    %Wl_bysim = cell(nsim,1);
    Wl_cv_bysim = cell(nsim,1);
    best_alpha_bysim = zeros(nsim,1);
    W_maskmat = load('W_mask.mat');
    W_mask = W_maskmat.W_mask;
    W_mask_active = W_mask(sum(W_mask~=0,2)~=0,sum(W_mask~=0,1)~=0);
    %criterion_cv_onfulldata_bysim = zeros(nsim,1);
    %criterion_old_onfulldata_bysim = zeros(nsim,1);
    for l = 1:nsim
        %Set the seed to be iterations number at the beginning of each iteration.
        rng(l);
        [~,~,~,W_init,~,X_input,S_input,~]=simdata(n,m,d,repprop,L_noise,X_noise,W_noise,W_mask_active);
        
        % Normalization of X
        X_input = normalize(X_input);
        
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
    [Wl_cv,~] = estW_smooth(X_input,X_input,S_input,W_init,best_alpha);
    %[Wl] = estW(X_input,S_input,W_init);

    fprintf('Completed for Simulation %i\n', l); 
    %X_bysim{l} = X;
    %X_input_bysim{l} = X_input;
    %W_bysim{l} = W;
    %W_init_bysim{l}  = W_init;
    %Wl_bysim{l} = Wl;
    Wl_cv_bysim{l} = Wl_cv;
    best_alpha_bysim(l) = best_alpha;
    %criterion_cv_onfulldata_bysim(l) = criterion_cv_onfulldata;
    %criterion_old_onfulldata_bysim(l) = criterion_old_onfulldata;
    disp('Results saved.');
    end

    fname = sprintf('results/cv/simres_all_n%d_m%d_d%d_repprop%.2f_Lnoise%.2f_Xnoise%0.2f_Wnoise%0.2f_nsim%d.mat', n, m, d, repprop, L_noise, X_noise, W_noise, nsim);
    save(fname,'Wl_cv_bysim','best_alpha_bysim','alphas','nfolds','-v7.3');
    disp('Success, All Results saved!');

end

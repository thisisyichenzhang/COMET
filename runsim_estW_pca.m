function runsim_estW_pca(n,m,d,repprop,L_noise,X_noise,W_noise,nsim)
    
    %X_bysim = cell(nsim,1);
    %X_input_bysim = cell(nsim,1);
    %W_bysim = cell(nsim,1);
    %W_init_bysim  = cell(nsim,1);
    Wl_pca_bysim = cell(nsim,1);
    L_pca_bysim = cell(nsim,1);
    W_maskmat = load('W_mask.mat');
    W_mask = W_maskmat.W_mask;
    W_mask_active = W_mask(sum(W_mask~=0,2)~=0,sum(W_mask~=0,1)~=0);
    %1832 198
    for l = 1:nsim
        %Set the seed to be iterations number at the beginning of each iteration.
        rng(l);
        [~,~,~,W_init,~,X_input,S_input,~]=simdata(n,m,d,repprop,L_noise,X_noise,W_noise,W_mask_active);


        [Wl_pca,L_pca] = estW_pca(X_input,S_input,W_init);

        fprintf('Completed for Simulation %i\n', l); 
    %X_bysim{l} = X;
    %X_input_bysim{l} = X_input;
    %W_bysim{l} = W;
    %W_init_bysim{l}  = W_init;
    Wl_pca_bysim{l} = Wl_pca;
    L_pca_bysim{l} = L_pca;
    disp('Results saved.');
    end

    fname = sprintf('results/estWpca/simres_estWpca_n%d_m%d_d%d_repprop%.2f_Lnoise%.2f_Xnoise%0.2f_Wnoise%0.2f_nsim%d.mat', n, m, d, repprop, L_noise, X_noise, W_noise, nsim);
    save(fname,'Wl_pca_bysim','L_pca_bysim','-v7.3');
    disp('Success, All Results saved!');

end
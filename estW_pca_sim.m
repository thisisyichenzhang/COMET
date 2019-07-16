function estW_pca_sim(n,m,d,repprop,L_noise,X_noise,W_noise,nsim,seed)
    %X_bysim = cell(nsim,1);
    %X_input_bysim = cell(nsim,1);
    %W_bysim = cell(nsim,1);
    %W_init_bysim  = cell(nsim,1);
    Wl_pca_bysim = cell(nsim,1);
    L_pca_bysim = cell(nsim,1);
    b = m/d;
    % Random seed.
    rng(seed);
    for l = 1:nsim
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

    fname = sprintf('results/estWpca/simres_estWpca_n%d_m%d_d%d_repprop%.2f_Lnoise%.2f_Xnoise%0.2f_Wnoise%0.2f_nsim%d_seed%d.mat', n, m, d, repprop, L_noise, X_noise, W_noise, nsim, seed);
    save(fname,'Wl_pca_bysim','L_pca_bysim','-v7.3');
    disp('Success, All Results saved!');

end
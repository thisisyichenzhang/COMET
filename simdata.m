clc; clear all;
for n = [100]
    paras.n = n;
    for m = [3000]
        paras.m = m;
        for d = [100]
            paras.d = d
            for S_noise = [0 1 10]
                paras.S_noise = S_noise;
                 for X_noise = [10]
                     paras.X_noise = X_noise;
                    for repprop = [0 0.25 0.8]
                        paras.repprop = repprop;
                         for nsim = [100]
                            paras.nsim = nsim;
                            [X_cell, W_cell, S_true, S_input, paras] = data(paras);
                            fname = sprintf('data/data_n%d_m%d_d%d_repprop%.1f_Snoise%.1f_Xnoise%d_nsim%d.mat', n, m, d, repprop, S_noise, X_noise, nsim);
                            save(fname);
                         end 
                    end
                 end
                 
            end
            
        end
        
    end
    
end


function [X_cell, W_cell, S_true, S_input, paras] = data(paras)
    n = paras.n;
    m = paras.m;
    d = paras.d;
    repprop = paras.repprop;
    X_noise = paras.X_noise;
    S_noise = paras.S_noise;
    nsim = paras.nsim;
    % Create empty cell array for X and W
    X_cell = cell(nsim,1);
    W_cell = cell(nsim,1);
    S_true_cell = cell(nsim,1);
    S_input_cell = cell(nsim,1);
    for j = 1:nsim
        % d/n  = 100/100; 500/100; 500/1000
        %m = 3000 and 15000 for d = 100 and 500 so that m/d ~= 30.
        %n = 1000; % Number of subjects
        %m = 100; % Number of Probes
        %d = 20;  % Number of genes
        %X_noise = 10 
        b = m/d; % Probes-Genes ratio
        X = zeros(n,m);

        W = zeros(m,d);
            for i = 1:d
                ind = (1+(i-1)*b):(i*b);
                if rand > repprop || i==1
                    sig = randn(n,1);
                end
                X(:,ind) = sig*ones(1,b) + X_noise*randn(n,b);
                W(ind,i) = rand(b,1);
            end
        X_cell{j} = X;
        W_cell{j} = W;
        L = X*W;
        
        S_true = corr(L);
        S_input = S_true + S_noise * randn(d,d); 
        S_true_cell{j} = S_true;
        S_input_cell{j} = S_input;
    end
    
end 
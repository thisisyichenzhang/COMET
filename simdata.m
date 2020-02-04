%% Function for generating data give W_mask (inferred from the real data)
% testing code
% rng(123);[X,W,L,W_init,L_input,X_input,S_input,paras] = simdata(500, 100, 20, 0.2, 0, 0, 0);
function[X,W,L,W_init,L_input,X_input,S_input,paras]=simdata(n,m,d,repprop,L_noise,X_noise,W_noise, W_mask)
        paras = struct; 
        paras.n = n; 
        paras.m = m; 
        paras.d = d;
        paras.repprop = repprop; 
        paras.L_noise = L_noise; 
        paras.X_noise = X_noise;
        paras.W_noise = W_noise;
        
        X = zeros(n,m);
        W = zeros(m,d);
        % Simulate the dataset
        for i = 1:d
            ind = logical(W_mask(:,i));
            number_of_active_probes= sum(ind);
            
            if rand > repprop || i==1
                sig = randn(n,1);
            end
            
                X(:,ind) = sig*ones(1,number_of_active_probes) + 1*randn(n,number_of_active_probes);
                W(ind,i) = rand(number_of_active_probes,1);
        
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
        X_input = X + X_noise .* randn(n,m);
            
end
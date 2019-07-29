clear;
n =100;
m = 4000;
d = 200;
X_noise = 50;
%L_noise = 0;
W_noise = 0;
repprop = 0.25;
b = m/d;

nsim = 1;
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
        %L_input = L + L_noise .* randn(n,d);
        
        X_input = X + X_noise .* randn(n,m);
        L_input = X_input * W_init;
        S_input = corr(L_input);
        
        %Normalization of X 
        %X_input = normalize(X_input);
        L_pca = estL_pca(X_input,W_init);
        L_init = X_input* W_init;
    end
    
mean(diag(corr(L,L_pca)))
mean(diag(corr(L,L_init)))
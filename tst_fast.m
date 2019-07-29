clear all;
n=100;
m=100;
d=20;
X_noise = 5;
L_noise= 0;
W_noise = 0;
nsim=1;
b = m/d;
repprop =0.5;

alpha=1;
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
        X_input = X + X_noise .* randn(n,m);
        [Wl,~] = estW_smooth_fast(X_input,X_input,S_input,W_init,alpha);
    end
    
Ll = X_input * Wl;  
Cl = corr(Ll);
L_init = X_input * W_init;
C_init = corr(L_init);
C = corr(L);
subplot(131);imagesc(C);subplot(132);imagesc(C_init);subplot(133);imagesc(Cl);
    
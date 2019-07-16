% Fixing L_pca
function [W,L_pca] = estW_pca(X_train,S,Wo)
    K = 1;% the number of iterations
    [n,~] = size(X_train);
    d = size(S,1);
    
% Initilization of W.    
%  W = randn(m,d);
%  W = Wo+randn(m,d);
%  W = double(Wo~=0);
%   W = double(Wo>0.25);
    W = Wo;
    
    % Compute L_PCA
    L_pca = zeros(n,d);
    for j=1:d
        fprintf('Computing PCA for gene %i\n', j);
        ind = W(:,j)~=0; 
        aaa = pca(X_train(:,ind)');
        L_pca(:,j) = aaa(:,1);
    end
    
    for k=1:K
        for i=1:d
            disp(i);
            options = optimoptions('fmincon','Display','iter');
            ind = W(:,i)~=0; 
            fun=@(w)(-1*(sum(([S(i,1:i-1) S(i,i+1:end)].*corr(X_train(:,ind)*w,[L_pca(:,1:i-1) L_pca(:,i+1:end)]))))/d);
            W(ind,i) = fmincon(fun,W(ind,i),[],[],[],[],[],[],[],options);
        end
    end
end


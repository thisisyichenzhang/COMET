% Fixing L_pca
function [L_pca] = estL_pca(X_train,Wo)

    [n,~] = size(X_train);
    d = size(Wo,2);
    
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
        [~,sc] = pca(X_train(:,ind));
        L_pca(:,j) = sc(:,1);
    end
    
end


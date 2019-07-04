% Function
function [W,criterion] = estW(X,S,Wo)
    K = 1; % Number of optimization iterations
%    [n,m] = size(X);
    d = size(S,1);
%   W = rand(m,d);
%   W = randn(m,d);
%   W = Wo+randn(m,d);
%   W = double(Wo~=0);
%   W = double(Wo>0.25);
    W = Wo;    
    L=X*W;
    for k=1:K
        for i=1:d
            disp(i);
            options = optimoptions('fmincon','Display','iter');
            ind = W(:,i)~=0;
            fun=@(w)(-1*sum(([S(i,1:i-1) S(i,i+1:end)].*corr(X(:,ind)*w,[L(:,1:i-1) L(:,i+1:end)])))/d);
            W(ind,i) = fmincon(fun,W(ind,i),[],[],[],[],[],[],[],options);
            L=X*W;
        end
    end
    criterion=mean(diag(transpose(W)*corr(X)*W)); 
end
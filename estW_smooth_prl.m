% use the correlation across probes to smooth the W. 
function [W,criterion] = estW_smooth_prl(X_train,X_test,S,Wo,alpha)
    K = 1;% the number of iterations
    %[n,m] = size(X);
    d = size(S,1);
    
% Initilization of W.    
%  W = randn(m,d);
%  W = Wo+randn(m,d);
%  W = double(Wo~=0);
%   W = double(Wo>0.25);
    W = Wo;
    L=X_train*W;
    for k=1:K
        for i=1:d
            disp(i);
            options = optimoptions('fmincon','Display','iter','MaxFunEvals',1e5);
            ind = W(:,i)~=0; 
            fun=@(w)(sum((corr(X_train(:,ind)*w,[L(:,1:i-1) L(:,i+1:end)]) - [S(i,1:i-1) S(i,i+1:end)]).^2)/(d-1) - alpha*sum(transpose(w)*corr(X_train(:,ind))*w)/(d-1));
            W(ind,i) = fmincon(fun,W(ind,i),[],[],[],[],[],[],[],options);
            % No update on L
            %L=X_train*W; 
        end
    end
    criterion=mean(diag(transpose(W)*corr(X_test)*W)); 
end


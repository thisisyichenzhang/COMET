% use the correlation across probes to smooth the W. 
function [W,criterion] = estW_smooth_new(X_train,X_test,S,Wo,alpha)
    K = 5;% the number of iterations
    %[n,m] = size(X);
    d = size(S,1);
    % Normilizing training X
    X_train = normalize(X_train);
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
            fun=@(w)(-w' * X_train(:,ind)' * L(:,i) - alpha * sum(L(:,i)' * X_train(:,ind) * w * [S(i,1:i-1) S(i,i+1:end)]));
            %fun=@(w)(w' * X_train(:,ind)' * X_train(:,ind) * w - alpha * sum([W(ind,1:i-1) W(ind,i+1:end)]' * X_train(:,ind)' * X_train(:,ind) * w * [S(i,1:i-1) S(i,i+1:end)]));
            %fun=@(w)(sum((corr(X_train(:,ind)*w,[L(:,1:i-1) L(:,i+1:end)]) - [S(i,1:i-1) S(i,i+1:end)]).^2)/(d-1) + alpha*sum(transpose(w)*corr(X_train(:,ind))*w)/(d-1));
            W(ind,i) = fmincon(fun,W(ind,i),[],[],[],[],[],[],[],options);
            % Update L
            L=X_train*W; 
        end
    end
    % Notice the nagtive sign, this is the likelihood to maximize.
    criterion= -mean(diag(W'*corr(X_test)*W)); 
end


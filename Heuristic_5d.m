% Code to predict the quasi-stationary infected fraction mu using the
% HEURISTIC 5d:
%
% \mu/n = (1/n) * \sum_{i=0}^{n} \hat{f}_i
%
% where \hat{f}_i is the solution to the system
%
% f_i = (\tau\sum_{j:j~i}\eta_{d_i,d_j}*f_j)/(1+\tau\sum_{j:j~i}\eta_{d_i,d_j}*f_j)
%
% with f_i and f_j the fraction of time node i and node j are
% respectively infected.
%
% \eta_{d_i,d_j} = P(Xj=1 | Xi=0)/P(Xj = 1)
%
% Input: 
% (i) n : The number of nodes. 
% (ii) tau : The infection rate
% (iii) EdgeList : Vector containing all neighbours of 1, then
%      all  neighbours of 2 and so on. Length: sum of degrees.
% (iv) Degr : Vector of length n giving degrees of all nodes.

% Output: 
% (i) Predicted infected fraction mu

function [Heuristic5d] = Heuristic_5d(n,tau,EdgeList,Degr)
    
    Adj_Matrix = zeros(n,n); %Adjacency matrix of the ER graph
    for i = 1:n
        d = EdgeList(EdgeList(:,1)==i,2);
        W = zeros(1,n);
        W(d)=1;
        Adj_Matrix(i,:) = W;
    end
    
    eta_matrix = zeros(n,n);
    [~,eta]=Heuristic_5b(n,tau,EdgeList,Degr);
    
    for i = 1:n
        for j = 1:n
            if Degr(i)*Degr(j) ~= 0
                eta_matrix(i,j) = eta(Degr(i),Degr(j));
            else
                eta_matrix(i,j) = 0;
            end

        end
    end
    
    % Define the initial guess for X
    Xvec = ones(n,1);
    
    % Define the tolerance and maximum number of iterations
    tol = 1e-6;
    maxiter = 10000;
    
    % convergence flag
    converged = false;
    
    % Iterate until convergence or maximum number of iterations
    for iter = 1:maxiter
        Xnew = zeros(n,1);
        for i = 1:n
            Xnew(i) = tau*sum(Xvec(Adj_Matrix(i,:)==1).*eta_matrix(i,[find(Adj_Matrix(i,:)==1)]))/(1 ...
                +tau*sum(Xvec(Adj_Matrix(i,:)==1).*eta_matrix(i,[find(Adj_Matrix(i,:)==1)])));
        end
        if norm(Xnew-Xvec) < tol
            converged = true;
            break;
        end
        Xvec = Xnew;
    end
    Heuristic5d = mean(Xvec);

end
% Code to predict the quasi-stationary infected fraction mu using the
% HEURISTICS :
%
% (i) \mu/n = \sum_{d=0}^{n-1} P(n,p,d) [(d\tau\mu)/n+d\tau\mu]
%
%     Here P(n,p,d) is the exact frequency #{i: d_i = d}/n, i.e., the
%     probability that a node i in the graph has degree d. The maximum degree 
%     here is ceil(n*p+5*sqrt(n*p*(1-p))). We will name this HEURISTIC 2b.

% (ii) \mu/n = \sum_{d=0}^{n-1} P(n,m/(nC2),d) [(d\tau\mu)/n+d\tau\mu]
%
%      Here P(n,p,d) is s the probability P(Bin(n-1,p)=d) or an approximation to it. 
%      p here is taken as the exact average degree m/(nC2). This heuristic 
%      uses the exact number of edges in the graph, which is equivalent to 
%      using the sum of the degrees. The maximum degree here is ceil(n*p+5*sqrt(n*p*(1-p))).
%      We will name this HEURISTIC 3.
%
% (iii) \mu/n = (1/n) * \sum_{i=0}^{n} \hat{f}_i
%
%       where \hat{f}_i is the solution to the system
%       f_i = (\tau\sum_{j:i ~ j}f_j)/(1+\tau\sum_{j:i ~ j}f_j)
%       with f_i and f_j the fraction of time node i and node j are
%       respectively infected. For the heuristic to be feasible, we assume
%       the full adjacency matrix of the Erdos Renyi (ER) graph is known and
%       there is an idependence between nodes. We will name this HEURISTIC
%       4.
%
% Input: 
% (i) n : The number of nodes. 
% (ii) tau : The infection rate
% (iii) EdgeList : Vector containing all neighbours of 1, then
%      all  neighbours of 2 and so on. Length: sum of degrees.
% (iv) Degr : Vector of length n giving degrees of all nodes.


% Output: 
% (i) Predicted infected fraction mu using HEURISTIC 2b
% (ii) Predicted infected fraction mu using HEURISTIC 3
% (iii) Predicted infected fraction mu using HEURISTIC 4

function [Heuristic_2b, Heuristic_3, Heuristic_4] = Heuristics2b_3_4(n,tau,EdgeList,Degr)
    % Here we calculate the predicted infected fraction mu using HEURISTIC 3. 
    % We use normal approximation to compute P(n,p,d) with p = 2m/(n(n-1))
    
    m = sum(Degr)/2;
    p = 2*m/(n*(n-1)); %exact average degree
    dmax = ceil(n*p+5*sqrt(n*p*(1-p))); % most likely larger than max degree
    step_3 = 0.1; % step-size 
    mu_3 = 0.5; % initial value of mu in HEURISTIC 3

    v = 0;
    for k_1 = 1:10
        if v < 1
            mu_3 = mu_3 - step_3;
        else
            mu_3 = mu_3 + step_3;
        end
        v = 0;
        for d = 1:min(n,dmax)
            % We calculate P(n,p,d) using the normal approximation 
            % P(|N((n-1)p,(n-1)p(1-p))-d|\leq 1/2)
            sigma = sqrt(p*(1-p)*(n-1));
            prob = normcdf(d+0.5,p*(n-1),sigma) - normcdf(d-0.5,p*(n-1),sigma);
            vd = prob*d*tau/(1+d*mu_3*tau);
            v = v+vd;
        end
        step_3 = (step_3)/2;
    end
    Heuristic_3 = mu_3;
    

    % Here we calculate the predicted infected fraction mu using HEURISTIC 2b. 
   
    [Freq, diff_degr] = groupcounts(Degr);
    probs = Freq/n; % probability of a node having a certain degree i.e. #{v:d(v)=d}/n
    step_2 = 0.1; % step-size 
    mu_2 = 0.5; % initial value of mu in HEURISTIC 2b
    
    y = 0;
    for k_2 = 1:10
        if y < 1
            mu_2 = mu_2 - step_2;
        else
            mu_2 = mu_2 + step_2;
        end
        y = 0;
        for l = 1:length(diff_degr)
            yl = sum(probs(l)*(diff_degr(l)*tau/(1+diff_degr(l)*tau*mu_2)));
            y = y+yl;
        end
        step_2 = (step_2)/2;
    end
    Heuristic_2b = mu_2;


    % Here we calculate the predicted infected fraction using HEURISTIC 4
    % We use the adjacency matrix to compute the infected fraction.

    Adj_Matrix = zeros(n,n); %Adjacency matrix of the ER graph
    for i = 1:n
        d = EdgeList(EdgeList(:,1)==i,2);
        W = zeros(1,n);
        W(d)=1;
        Adj_Matrix(i,:) = W;
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
            Xnew(i) = tau*sum(Xvec(Adj_Matrix(i,:)==1))/(1+tau*sum(Xvec(Adj_Matrix(i,:)==1)));
        end
        if norm(Xnew-Xvec) < tol
            converged = true;
            break;
        end
        Xvec = Xnew;
    end
    Heuristic_4 = mean(Xvec);
end
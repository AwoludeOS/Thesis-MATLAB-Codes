% Code to predict the quasi-stationary infected fraction mu using the
% HEURISTIC 2a:
%
% \mu/n = \sum_{d=0}^{n-1} P(n,p,d) [(d\tau\mu)/n+d\tau\mu]
%
% Here P(n,p,d) is the degree distribution of each node i.e., it is the
% probability P(Bin(n-1,p)=d). We can also interpret this as the probability 
% that a node i in the graph has degree d (d here is possible degree of 
% a node in the graph). The maximum degree here is ceil(n*p+5*sqrt(n*p*(1-p))).
%
% Input: 
% (i) n : The number of nodes. 
% (ii) p : Edge probability.
% (iii) tau : The infection rate

% Output: 
% (i) Predicted infected fraction mu

function [Hueristic_2a] = Heuristic2a_bin(n,p,tau)
    dmax = ceil(n*p+5*sqrt(n*p*(1-p))); % most likely larger than max degree
    step = 0.1; % step-size 
    mu = 0.5; % initial value of mu
    
    % Here we calculate in infected fraction using the heuristic stated
    % above
    u = 0;
    for k = 1:10
        if u < 1
            mu = mu - step;
        else
            mu = mu + step;
        end
        u = 0;
        for d = 1:min(n,dmax)
            % We calculate P(n,p,d) using the probability P(Bin(n-1,p)=d)
            prob = binopdf(d, n-1, p);
            ud = prob*d*tau/(1+d*mu*tau);
            u = u+ud;
        end
        step = (step)/2;
    end
    Hueristic_2a = mu;
end
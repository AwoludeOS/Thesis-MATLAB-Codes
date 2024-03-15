% Code to predict the quasi-stationary infected fraction mu using the
% HEURISTIC 5c:
%
% \mu/n = \sum_{d_i=1}^{n-1}\sum_{d_j=1}^{n-1} P(n,p,d_i)*P(n,p,d_j)*[x_3(d_i,d_j)+x_4(d_i,d_j)]
%
% Here P(n,p,d)=P^{\geq 1}(n,p,d) and P(n,p,e) are the probabilities that 
% node i and node j have degrees d and e respectively. That is 
%
% P(n,p,d)*P(n,p,e) = [1/(#{i:d_i\geq 1})]*[\sum_{i=1}^{n}\sum_{j~i} (1{d_i=d, d_j=e})/d)]
%
% The maximum degree here is ceil(n*p+5*sqrt(n*p*(1-p))).
%
% Input: 
% (i) n : The number of nodes. 
% (ii) tau : The infection rate
% (iii) EdgeList : Vector containing all neighbours of 1, then
%      all  neighbours of 2 and so on. Length: sum of degrees.
% (iv) Degr : Vector of length n giving degrees of all nodes.

% Output: 
% (i) Predicted infected fraction mu

function [heuristic5c] = Heuristic_5c(n,tau,EdgeList,Degr)
    
    dmax = max(Degr);

    prob = zeros(dmax,dmax);

    for i = 1:length(EdgeList)
        d = Degr(EdgeList(i,1));
        e = Degr(EdgeList(i,2));
        prob(d,e) = prob(d,e)+(1/d);
    end
    prob = prob/sum(Degr>=1);

    
    % observe that sum(pj.*diff_degr) is approximately np+1 which 
    % is the average degree of node j and sum(pi.*diff_degr) is 
    % approximately np which is the average degree of node i
    % 'Maximal' degree: expectation + 5 standard deviations

    epsilon = 10^(-5);
    
    museta = 0;
    musco = 0.5;
    
    
    while abs(museta-musco) > epsilon
        museta = musco;
    
        mu = zeros(dmax,dmax);
        mus = zeros(dmax,dmax);
        mucond = zeros(dmax,dmax);
        muscond = zeros(dmax,dmax);
        eta = zeros(dmax,dmax);
        
    
        for di = 1:dmax
            for dj = 1:dmax
                % Consider a Markov chain on 4 states:
                % (Xi,Xj) = (0,0), (0,1), (1,0), (1,1).
                % Generator matrix transition rates:
                Q = [0 (dj-1)*museta*tau (di-1)*museta*tau 0;
                    1 0 0 tau*(di-1)*museta+tau;
                    1 0 0 (dj-1)*museta*tau+tau;
                    0 1 1 0];
                % Generator matrix diagonal:
                Q = Q-diag(sum(Q,2));
                % solve for stationary distribution:
                q = null(Q');
                q = q/sum(q);
                q1 = q(1);q2 = q(2);q3 = q(3);q4 = q(4);
    
                % probabilities of infection for i and j:
                mu(di,dj) = q3+q4;
                mus(di,dj) = q2+q4;
                % conditional probabilities P(Xi=1 | Xj=0) and P(Xj=1 | Xi=0):
                mucond(di,dj) = q3/(1-mus(di,dj));
                muscond(di,dj) = q2/(1-mu(di,dj));

                % approximation for eta:
                
                if di*dj ~=1
                    eta(di,dj) = muscond(di,dj)/mus(di,dj); % this one seems most sensible to use
                else
                    eta(di,dj) = 0.000001;
                end
    
           
            end
            
        end
        musco = sum(sum(prob.*muscond));
    end
    
    heuristic5c = (sum(Degr>=1)*sum(sum(prob.*mu)))/n;
    
end
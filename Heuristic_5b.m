% Code to predict the quasi-stationary infected fraction mu using the
% HEURISTIC 5b :
%
% \mu/n = \sum_{d_i=1}^{n-1}\sum_{d_j=1}^{n-1} P(n,p,d_i)*P(n,p,d_j)*[x_3(d_i,d_j)+x_4(d_i,d_j)]
%
% Here P(n,p,d_i) and P(n,p,d_j) are the probabilities that node i and node
% j have degrees d_i and d_j respectively. That is 
% P(n,p,d_i) = (1/n)*\sum_{i=1}^{n} 1{d_i=d},
% P% Code to predict the quasi-stationary infected fraction mu using the
% heuristic :
%
% \mu/n = \sum_{d_i=1}^{n-1}\sum_{d_j=1}^{n-1} P(n,p,d_i)*P(n,p,d_j)*[x_3(d_i,d_j)+x_4(d_i,d_j)]
%
% Here P(n,p,d_i) and P(n,p,d_j) are the probabilities that node i and node
% j have degrees d_i and d_j respectively. That is 
%
% P(n,p,d_i) = P^{\geq 1}(n,p,d_i) = (\sum_{i=1}^{n} 1{d_i=d})/(\sum_{i=1}^{n} 1{d_i\geq 1}) and 
%
% P(n,p,d_j) = [1/(#{i:d_i\geq 1})]*[\sum_{i=1, d_i\geq 1}^{n}(1/d_i) * (\sum_{j~i} 1{d_i=d})].
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
% (i) Predicted infected fraction mu(n,p,d_i) = (\sum_{i=1}^{n} 1{d_i=d})/(\sum_{i=1}^{n} 1{d_i\geq 1}), and 
% P(n,p,d_j) = \sum_{i=1, d_i\geq 1}^{n}(1/d_i) * (\sum_{j~i} 1{d_i=d}).
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

function [heuristic5b,eta] = Heuristic_5b(n,tau,EdgeList,Degr)
    
    % We calculate the probability P(n,p,d_i)
    pi = zeros(max(Degr),1); 
    for d = 1:max(Degr)
        pi(d) = sum(Degr==d)/sum(Degr>=1);
    end

    % We calculate the probability P(n,p,d_j)
    Degr_Pairs = Degr(EdgeList);
    Frac_Degr_Pairs = 1./Degr_Pairs(:,1);
    pj = zeros(max(Degr),1); 
    
    % We calculate the probability P(n,p,d_i)*P(n,p,d_j)
    prob = zeros(max(Degr),max(Degr));
    for i = 1:max(Degr)
        for j = 1:max(Degr)
            pj(j) = sum(Frac_Degr_Pairs(Degr_Pairs(:,2)==j))/(sum(Degr>=1));
            prob(i,j) = pi(i)*pj(j);  
        end
    end

    % Observe that sum(pj.*diff_degr) is approximately np+1 which 
    % is the average degree of node j and sum(pi.*diff_degr) is 
    % approximately np which is the average degree of node i
    % 'Maximal' degree: expectation + 5 standard deviations

    
    epsilon = 10^(-5);
    
    museta = 0;
    musco = 0.5;
    
    dmax = max(Degr);
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
                % eta is the estimate of the ratio of the probabilities
                % P(Xj=1 | Xi=0) and P(Xj=1)
                if di*dj ~=1
                    eta(di,dj) = muscond(di,dj)/mus(di,dj); 
                else
                    eta(di,dj) = 0.000001;
                end
    
           
            end
            
        end
        musco = sum(sum(prob.*muscond));
    end
    
    heuristic5b = (sum(Degr>=1)*sum(sum(prob.*mu)))/n;
    
end
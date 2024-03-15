% Code to predict the quasi-stationary infected fraction mu using the
% HEURISTIC 5 :
%
% \mu/n = \sum_{d_i=1}^{n-1}\sum_{d_j=1}^{n-1} P(n,p,d_i)*P(n,p,d_j)*[x_3(d_i,d_j)+x_4(d_i,d_j)]
%
% Here P(n,p,d_i) and P(n,p,d_j) are the probabilities that node i and node
% j have degrees d_i and d_j respectively. That is P(n,p,d_i) = P(B = d_i | B \geq 1)
% and P(n,p,d_j) = P(\tilde{B} = d_j) with B ~ Bin(n-1,p) and \tilde{B} ~ Bin(n-2,p).
% The maximum degree here is ceil(n*p+5*sqrt(n*p*(1-p))).
%
% Input: 
% (i) n : The number of nodes. 
% (ii) p : Edge probability.
% (iii) tau : The infection rate

% Output: 
% (i) Predicted infected fraction mu

function Heur5_mu = Heuristic5(n,p,tau)

    % 'Maximal' degree: expectation + 5 standard deviations
    dmax = ceil(n*p+5*sqrt(n*p*(1-p)));

    kansi0 = (1-p)^(n-1);
    pi = zeros(dmax,1); % This is P(n,p,d_i)
    pj = zeros(dmax,1); % This is P(n,p,d_j)
    prob = zeros(dmax,dmax); %This is P(n,p,d_i)*P(n,p,d_j)
    
    % Here we calculate the probability P(n,p,d_i)*P(n,p,d_j)
    for di = 1:dmax
        for dj = 1:dmax
            pi(di) = nchoosek(n-1,di)*p^di*(1-p)^(n-1-di)/(1-kansi0); % P(deg(i) = di)
            pj(dj) = nchoosek(n-2,dj-1)*p^(dj-1)*(1-p)^(n-2-(dj-1)); % P(deg(j) = dj)
    
            prob(di,dj) = pi(di)*pj(dj);  
        end
    end

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

                % conditional probabilities P(Xi=1 | Xj=0) and P(Xi=1 | Xj=0):
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
    
    Heur5_mu = (1-kansi0)*sum(sum(prob.*mu));
    
end
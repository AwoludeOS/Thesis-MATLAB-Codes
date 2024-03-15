% Code to simulate the contact process on the Erdos Renyi (ER) Graph.
% We first generate an ER graph with parameters n and p. Then run the contact
% process with infection rate tau and healing rate 1 during an interval of
% length time. At t=0, nodes are healthy with probability 1/(n*p*tau) and
% infected otherwise. This is the theoretical healthy fraction in the
% limit. If n*p*tau = 2, the code is efficient. 

% Input: 
% (i) n : The number of nodes. 
% (ii) p : Edge probability.
% (iii) tau : The infection rate
% (iv) time : Time for the contact process run on the ER graph

% Output: 
% (i) Time_Per_State : The time spent in each state (state = number of
% infected).
% (ii) EdgeList : Vector containing all neighbours of 1, then
%      all  neighbours of 2 and so on. Length: sum of degrees.
% (iii) Degr : Vector of length n giving degrees of all nodes.


function [Time_Per_State,EdgeList,Degr] = CP_on_Erdos_Renyi_V1(n,p,tau,time)

    [Nei_Mat,EdgeList,Degr] = Erdos_Renyi(n,p); % We generate the ER graph using
    % parameters n and p

    no_edges = size(EdgeList,1);% The number of edges in the ER graph
    
    X=rand(n,1)<(n-1/(p*tau))/n;% List of infectuous states at time 0
    
    Heal_rate = sum(X);% Healing rate in the process
    
    no_edges_active = sum(X(EdgeList(:,1)).*(1-X(EdgeList(:,2)))); % We keep track the number 
    % of active edges as the process evolve.

    Inf_rate = tau*no_edges_active; % Infection rate in the process
    
    T = 0; % We set a counter for time for the process to run
    
    Time_Per_State = zeros(n,1); % A vector that keeps track of the time spent per state.

    % Each entry in the vector "nodes_inftime_list"
    % represents the total time each node remained infected.
    
    while (T < time) && (Heal_rate > 0)
        
        dt=exprnd(1/(Inf_rate + Heal_rate)); % Exponential time taken for 
        % an event to occur

        T=T+dt;

        Time_Per_State(Heal_rate)=Time_Per_State(Heal_rate)+dt;
        
    
        Inf_Prob=Inf_rate/(Inf_rate+Heal_rate); % Infection probability
    
        if rand < Inf_Prob
            % We will get an infection
            active = 0;
            while active == 0
                index = ceil(rand*no_edges); % Pick a random edge
                node1 = EdgeList(index,1);
                node2 = EdgeList(index,2);
                % Check if the edge (node1, node2) is active
                active = X(node1)*(1-X(node2));
            end
           
            % node1 is going to infect node2
            X(node2)=1;
            Heal_rate=Heal_rate+1;
            % Update no_edges_active
            for i = 1:Degr(node2)% Check the neighbours of node2
                status = X(Nei_Mat(node2,i)); % Status of neighbour
                no_edges_active = no_edges_active +(1-2*status);
            end
            
            Inf_rate = tau*no_edges_active; % We update the infection rate 
            % based on the number of active edges
    
        else
            % We will get an healing
            infected = 0;
            while infected == 0
                node = ceil(rand*n); % We pick a random node
                infected = X(node);
            end
    
            X(node) = 0;
            Heal_rate=Heal_rate-1;
            % We update no_edges_active
            for i = 1:Degr(node)% Check the neighbours of node
                status = X(Nei_Mat(node,i)); % Status of neighbour
                no_edges_active = no_edges_active + (2*status-1);
            end
           
            Inf_rate = tau*no_edges_active;% We update the infection rate 
            % based on the number of active edges
    
        end
      
    end
end


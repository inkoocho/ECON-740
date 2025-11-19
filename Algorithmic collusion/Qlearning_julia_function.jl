##
##    Written by Jung Jae Kim on January 10, 2025. Translated from Qlearning_duopoly_dist_GR.m file
##

function Qlearning_julia_function(n_periods,n_prices,kbar,alpha,delta,epsilon)

    # Set up the game parameters
    #n_periods = 10000;
    n_players = 2;
    #n_prices = 15;
    n_init = 2;
    #kbar = 5; # exploration neighborhood parameter 1 ~ 14

    #alpha = 0.05;                          # Learning rate
    #beta = 2*10^(-5);                      # Epsilon-greedy parameter
    #delta = 0.95;                          # Discount factor
    eta = 0.1;                              # Prices range parameter

    # Set up demand function parameters
    # actual demand q_i=Ai-Bi p_i +Ci p_j   B>C
    A=1;
    B=1;
    C=0.7;

    A1=A;
    A2=A;
    B1=B;
    B2=B;
    C1=C;
    C2=C;

    # Set up prices range
    p_nash = (2*A*B+ C*A)/(4*B^2 - C^2);
    p_collude = (2*A*B+2*A*C) / (4*(B^2-C^2));
    price_range = LinRange(p_nash - eta*(p_collude-p_nash),p_collude + eta*(p_collude-p_nash),n_prices);

    # Initialize Q-matrix for each player 
    Q_0 = zeros(n_prices, n_prices, n_prices, n_players, n_init);

    # Q_0 for Player 1 (Collusive highest initial Q-matrix)
    for m1 = 1:n_prices
        p = [price_range[m1], price_range[m1]];
        q = [max(A1-B1*p[1]+C1*p[2],0), max(A2-B2*p[2]+C1*p[1],0)];
        payoff = p.*q;
        Q_0[:,:,m1,1,2] .= payoff[1] / (1-delta);
    end

    # Q_0 for Player 2 (Collusive highest initial Q-matrix)
    for m2 = 1:n_prices
        p = [price_range[m2], price_range[m2]];
        q = [max(A1-B1*p[1]+C1*p[2],0), max(A2-B2*p[2]+C1*p[1],0)];
        payoff = p.*q;
        Q_0[:,:,m2,2,2] .= payoff[1] / (1-delta);
    end

    # Nash highest initial Q-matrix
    Q_0[:,:,:,1,1] = reverse(Q_0[:,:,:,1,2],dims=3);
    Q_0[:,:,:,2,1] = reverse(Q_0[:,:,:,2,2],dims=3);

    prices = zeros(n_players, n_periods, n_init);
    optimal_prices = zeros(n_players, n_periods, n_init);
    quantity = zeros(n_players, n_periods, n_init);
    payoff = zeros(n_players, n_periods, n_init);
    joint_dist = zeros(n_prices,n_prices, n_init);
    Q_vec = zeros(n_prices,n_players, n_init);

    # Initial state
    Q = Q_0;
    state1 = [rand(1:n_prices),rand(1:n_prices)];
    prices[:,1,1] = [price_range[state1[1]],price_range[state1[2]]];
    optimal_prices[:,1,1] = prices[:,1,1];
    joint_dist[state1[1],state1[2],1] = 1;

    state2 = [rand(1:n_prices),rand(1:n_prices)];
    prices[:,1,2] = [price_range[state2[1]],price_range[state2[2]]];
    optimal_prices[:,1,2] = prices[:,1,2];
    joint_dist[state2[1],state2[2],2] = 1;
    
    converge_period = 1;

    for t = 2:n_periods
        #epsilon = 0.1;
        #epsilon = exp(-beta*t);  % time decreasing epsilon-greedy

        # Calculate Q-vector: Q-matrix * conditional distribution
        for s = 1:n_prices
            Q_vec[s,1,1] = (reshape(Q[state1[1],:,s,1,1],(1,n_prices)) * 
                            reshape((joint_dist[state1[1],:,1]/sum(joint_dist[state1[1],:,1])),(n_prices,1)))[];
            Q_vec[s,2,1] = (reshape(Q[:,state1[2],s,2,1],(1,n_prices)) *
                            reshape((joint_dist[:,state1[2],1]/sum(joint_dist[:,state1[2],1])),(n_prices,1)))[];

            Q_vec[s,1,2] = (reshape(Q[state2[1],:,s,1,2],(1,n_prices)) *
                            reshape((joint_dist[state2[1],:,2]/sum(joint_dist[state2[1],:,2])),(n_prices,1)))[];
            Q_vec[s,2,2] = (reshape(Q[:,state2[2],s,2,2],(1,n_prices)) *
                            reshape((joint_dist[:,state2[2],2]/sum(joint_dist[:,state2[2],2])),(n_prices,1)))[];
        end
    
        # Loop over players
        for i = 1:n_players
            # Epsilon-greedy exploration
            if rand() < epsilon
                # Choose a random price
                # prices(i,t) = price_range(randi(n_prices));
    
                # Choose a random price from the previous price state
                exploration_range1 = price_range[max(1,state1[i][1]-kbar):min(15,state1[i][1]+kbar)];
                prices[i,t,1] = exploration_range1[rand(1:length(exploration_range1))];
            else          
                # Choose the best price based on Q-vector
                idx1 = findmax(Q_vec[:,i,1])[2];
                optimal_prices[i,t,1] = price_range[idx1];
                prices[i,t,1] = price_range[idx1];
            end
        end

        for i = 1:n_players
            # Epsilon-greedy exploration
            if rand() < epsilon
                # Choose a random price
                # prices(i,t) = price_range(randi(n_prices));

                # Choose a random price from the previous price state
                exploration_range2 = price_range[max(1,state2[i][1]-kbar):min(15,state2[i][1]+kbar)];
                prices[i,t,2] = exploration_range2[rand(1:length(exploration_range2))];
            else
                # Choose the best price based on Q-vector
                idx2 = findmax(Q_vec[:,i,2])[2];
                optimal_prices[i,t,2] = price_range[idx2];
                prices[i,t,2] = price_range[idx2];
            end
        end

        #Compute quantities and payoffs for each player
        quantity[1,t,1] = max(A1-B1*prices[1,t,1]+C1*prices[2,t,1],0);
        quantity[2,t,1] = max(A2-B2*prices[2,t,1]+C2*prices[1,t,1],0);

        payoff[1,t,1] = prices[1,t,1]*quantity[1,t,1];
        payoff[2,t,1] = prices[2,t,1]*quantity[2,t,1];

        quantity[1,t,2] = max(A1-B1*prices[1,t,2]+C1*prices[2,t,2],0);
        quantity[2,t,2] = max(A2-B2*prices[2,t,2]+C2*prices[1,t,2],0);

        payoff[1,t,2] = prices[1,t,2]*quantity[1,t,2];
        payoff[2,t,2] = prices[2,t,2]*quantity[2,t,2];

        # Compute indices of chosen prices
        tolerance = 1e-5  # Adjust the tolerance value as needed
        action1 = [findall(x-> isapprox(x, prices[1,t,1], atol=tolerance), price_range),
                   findall(x-> isapprox(x, prices[2,t,1], atol=tolerance), price_range)];

        action2 = [findall(x-> isapprox(x, prices[1,t,2], atol=tolerance), price_range),
                   findall(x-> isapprox(x, prices[2,t,2], atol=tolerance), price_range)];

        #Update joint distribution and Q-value for each player
        p_mat1 = zeros(n_prices,n_prices);
        p_mat1[action1[1],action1[2]] .= 1;
        joint_dist[:,:,1] = (1-alpha)*joint_dist[:,:,1] + alpha*p_mat1;
        # joint_dist = ((t-1)/t)*joint_dist + (1/t)*p_mat;

        p_mat2 = zeros(n_prices,n_prices);
        p_mat2[action2[1],action2[2]] .= 1;
        joint_dist[:,:,2] = (1-alpha)*joint_dist[:,:,2] + alpha*p_mat2;  

        for i = 1:n_players
            Q[state1[1],state1[2], action1[i], i,1] = (1-alpha) * Q[state1[1],state1[2], action1[i], i,1] .+ 
                alpha * (payoff[i,t,1] + delta * maximum(Q[action1[1],action1[2], :, i,1]));

            Q[state2[1],state2[2], action2[i], i,2] = (1-alpha) * Q[state2[1],state2[2], action2[i], i,2] .+
                alpha * (payoff[i,t,2] + delta * maximum(Q[action2[1],action2[2], :, i,2]));
        end

        # Update state
        state1 = action1;
        state2 = action2;

        # Gelman-Rubin Statistic
        if mod(t,1000000)==0
            n = 1000000;
            W = zeros(2,2);
            B = zeros(2,2);

            avg_1 = mean(prices[:,t-999999:t,1],2);
            avg_2 = mean(prices[:,t-999999:t,2],2);
            avg_all = mean([avg_1,avg_2],2);

            W[:,:] = (1/(2*(n-1)))*((prices[:,t-999999:t,1]-avg_1)*(prices[:,t-999999:t,1]-avg_1)' + 
                        (prices[:,t-999999:t,2]-avg_2)*(prices[:,t-999999:t,2]-avg_2)');
            B[:,:] = (avg_1-avg_all)*(avg_1-avg_all)' + (avg_2-avg_all)*(avg_2-avg_all)';
            GR = (n-1)/n + (3/2)*eigs(inv(W[:,:])*B[:,:],1);

            if GR < 1.005
                break
            end
        end
        
        converge_period = converge_period + 1
    end    
    
    return prices, converge_period
end
function output=Qlearning_duopoly_dist_fun()
% Set up the game parameters
n_periods = 20000000;
n_players = 2;
n_prices = 15;
kbar = 5; % exploration neighborhood parameter 1 ~ 14

alpha = 0.04;                           % Learning rate
%alpha = 0.15;  
%beta = 2*10^(-5);                       % Epsilon-greedy parameter
delta = 0.95;                           % Discount factor
eta = 0.1;                              % Prices range parameter

% Set up demand function parameters
%actual demand q_i=Ai-Bi p_i +Ci p_j   B>C
A=1;
B=1;
C=0.7;

A1=A;
A2=A;
B1=B;
B2=B;
C1=C;
C2=C;

% Set up prices range
p_nash = (2*A*B+ C*A)/(4*B^2 - C^2);
p_collude = (2*A*B+2*A*C) / (4*(B^2-C^2));
price_range = linspace(p_nash - eta*(p_collude-p_nash),p_collude + eta*(p_collude-p_nash),n_prices);

% Initialize Q-matrix for each player (equation (8) in the page 3275)
Q_0 = zeros(n_prices, n_prices, n_prices, n_players);
% Q_0 for Player 1
for m1 = 1:n_prices
    payoff = [0, 0];
    for m2 = 1:n_prices
        p = [price_range(m1), price_range(m2)];
        q = [max(A1-B1*p(1)+C1*p(2),0), max(A2-B2*p(2)+C1*p(1),0)];
        payoff = payoff + p.*q;
    end
    Q_0(:,:,m1,1) = payoff(1) / ((1-delta)*n_prices);
end

% Q_0 for Player 2
for m2 = 1:n_prices
    payoff = [0, 0];
    for m1 = 1:n_prices
        p = [price_range(m1), price_range(m2)];
        q = [max(A1-B1*p(1)+C1*p(2),0), max(A2-B2*p(2)+C1*p(1),0)];
        payoff = payoff + p.*q;
    end
    Q_0(:,:,m2,2) = payoff(2) / ((1-delta)*n_prices);
end

% % Q_0 for Player 1
% for m1 = 1:n_prices
%     p = [price_range(m1), price_range(m1)];
%     q = [max(A1-B1*p(1)+C1*p(2),0), max(A2-B2*p(2)+C1*p(1),0)];
%     payoff = p.*q;
%     Q_0(:,:,m1,1) = payoff(1) / (1-delta);
% end
% 
% % Q_0 for Player 2
% for m2 = 1:n_prices
%     p = [price_range(m2), price_range(m2)];
%     q = [max(A1-B1*p(1)+C1*p(2),0), max(A2-B2*p(2)+C1*p(1),0)];
%     payoff = p.*q;
%     Q_0(:,:,m2,2) = payoff(1) / (1-delta);
% end

prices = zeros(n_players, n_periods);
optimal_prices = zeros(n_players, n_periods);
quantity = zeros(n_players, n_periods);
payoff = zeros(n_players, n_periods);
joint_dist = zeros(n_prices,n_prices);
Q_vec = zeros(n_prices,n_players);

% Initial state
Q = Q_0;
state = [randi(n_prices),randi(n_prices)];
prices(:,1) = [price_range(state(1)),price_range(state(2))];
optimal_prices(:,1) = prices(:,1);
joint_dist(state(1),state(2)) = 1;

% Loop over periods
for t = 2:n_periods
    epsilon = 0.1;
    %epsilon = exp(-beta*t);  % time decreasing epsilon-greedy

    % Calculate Q-vector: Q-matrix * conditional distribution
    for s = 1:n_prices
        Q_vec(s,1) = Q(state(1),:,s,1) * (joint_dist(state(1),:)/sum(joint_dist(state(1),:)))';
        Q_vec(s,2) = Q(:,state(2),s,2)' * (joint_dist(:,state(2))/sum(joint_dist(:,state(2))));
    end

    % Loop over players
    for i = 1:n_players
        % Epsilon-greedy exploration
        if rand < epsilon
            % % Choose a random price
            % prices(i,t) = price_range(randi(n_prices));

            % Choose a random price from the previous price state
            exploration_range = price_range(max(1,state(i)-kbar):min(15,state(i)+kbar));
            prices(i,t) = exploration_range(randi(length(exploration_range)));
        else          
            % Choose the best price based on Q-vector
            [~,idx] = max(Q_vec(:,i));
            optimal_prices(i,t) = price_range(idx);
            prices(i,t) = price_range(idx);
        end
    end

    % Compute quantities and payoffs for each player
    quantity(1,t) = max(A1-B1*prices(1,t)+C1*prices(2,t),0);
    quantity(2,t) = max(A2-B2*prices(2,t)+C2*prices(1,t),0);

    payoff(1,t) = prices(1,t)*quantity(1,t);
    payoff(2,t) = prices(2,t)*quantity(2,t);

    % Compute indices of chosen prices
    action = [find(price_range == prices(1,t)),find(price_range == prices(2,t))];

    % Update joint distribution and Q-value for each player
    p_mat = zeros(n_prices,n_prices);
    p_mat(action(1),action(2)) = 1;
    joint_dist = (1-alpha)*joint_dist + alpha*p_mat;
    % joint_dist = ((t-1)/t)*joint_dist + (1/t)*p_mat;

    for i = 1:n_players
        Q(state(1),state(2), action(i), i) = (1-alpha) * Q(state(1),state(2), action(i), i) + ...
            alpha * (payoff(i,t) + delta * max(Q(action(1),action(2), :, i)));
    end

    % Update state
    state = action;
end

output = prices(:,end-499999:end);
%output = [prices; optimal_prices];
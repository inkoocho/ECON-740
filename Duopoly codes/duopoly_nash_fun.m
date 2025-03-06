function [p_nash1,p_nash2,br, profit]=duopoly_nash_fun(N,lambda,wn)
% Model parameters
%N = 10000;

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

% NE
%p_nash = A / (2*B-C);
%p_collude=A /(2*(B-C));

p_nash1=(2*A1*B2+ C1*A2)/(4*B1*B2 - C1*C2);
p_nash2=(2*A2*B1+ C2*A1)/(4*B1*B2 - C1*C2);

q_nash1=A1-B1*p_nash1+C1*p_nash2;
q_nash2=A2-B2*p_nash2+C2*p_nash1;

p_collude1= (2*A1*B2+A2*(C1+C2)) / (4*B1*B2-(C1+C2)*(C1+C2));
p_collude2= (2*A2*B1+A1*(C1+C2)) / (4*B1*B2-(C1+C2)*(C1+C2));

q_collude1=A1-B1*p_collude1+C1*p_collude2;
q_collude2=A2-B2*p_collude2+C2*p_collude1;

p_leader= (2*A1*B2+ C1*A2)/(4*B1*B2 - 2*C1*C2);
p_follower = (A2 + C2*p_leader)/(2*B2);

q_leader=A1-B1*p_leader+C1*p_follower;
q_follower=A1-B1*p_follower+C1*p_leader;


% set the matrix size
p = zeros(2,N);

%alpha1 = zeros(2,N);
alpha1 = zeros(1,N);
%alpha2 = zeros(2,N);
alpha2 = zeros(1,N);

beta1=zeros(3,N);
beta2=zeros(3,N);

btemp1 =zeros(2);
btemp2 =zeros(2);

br=zeros(2,N);

pavg=zeros(2,N);

profit=zeros(2,N);

% sigmawn=0.065;
% wn = sigmawn*randn(2,N);

%R1=eye(2,2);
%R2=eye(2,2);

RR1=eye(3,3);
RR2=eye(3,3);

% initialize

%alpha1(:,1) = [p_nash1;0];
alpha1(:,1) = p_nash1;
%alpha2(:,1) = [p_nash2;0];
alpha2(:,1) = p_nash2;

beta1(:,1)=[A1; B1; C1];
beta2(:,1)=[A2; B2; C2];

%sig = [0.1, 0.1];
%lambda = 0.0025;

p(:,1) = [p_nash1; p_nash2];
pavg(:,1)=p(:,1);

for t = 1:N-1

    br(1,t)=(beta1(1,t)+ beta1(3,t)*alpha1(1,t)) / (2*(beta1(2,t)));
    br(2,t)=(beta2(1,t)+ beta2(3,t)*alpha2(1,t)) / (2*(beta2(2,t)));

    p(1,t) = br(1,t) +wn(1,t);
    p(2,t) = br(2,t) +wn(2,t);

    q1 = max(0,A1 - B1*p(1,t) + C1*p(2,t));
    q2 = max(0,A2 - B2*p(2,t) + C2*p(1,t));

    profit(:,t) = [p(1,t)*q1; p(2,t)*q2];

    pavg(:,t+1)=pavg(:,t) +lambda*( p(:,t) -pavg(:,t));

    % R1=[ 1, br(1,t); br(1,t), br(1,t)*br(1,t)]+diag([0, wn(1,t)*wn(1,t)]);
    alpha1(1,t+1) = alpha1(1,t) + lambda * (p(2,t) - alpha1(1,t));
    % R1 = R1 +lambda*( [1 , p(1,t); p(1,t), p(1,t)*p(1,t)] -R1);

    %upper= max((p_collude2 - p_nash2)/(p_collude1 - p_nash1) ,(p_collude1-p_nash1) / (p_collude2 - p_nash2));

    %if alpha1(2,t+1) > upper
    %    alpha1(2,t+1)= upper ;
    %    alpha1(1,t+1)=pavg(2,t+1)-pavg(1,t+1);
    %elseif alpha1(2,t+1) <0
    %    alpha1(2,t+1)=0;
    %    alpha1(1,t+1)=pavg(2,t+1);
    %end


    %R2=[ 1, br(2,t); br(2,t), br(2,t)*br(2,t)]+diag([0, wn(2,t)*wn(2,t)]);
    alpha2(1,t+1) = alpha2(1,t) + lambda * (p(1,t) - alpha2(1,t));
    %R2 = R2 +lambda*( [1 , p(2,t); p(2,t), p(2,t)*p(2,t)] -R2);

    %if alpha2(2,t+1) > upper
    %    alpha2(2,t+1)=upper;
    %    alpha2(1,t+1)=pavg(1,t+1)-pavg(2,t+1);
        %elseif alpha2(2,t+1) <0
        %    alpha2(2,t+1)=0;
        %    alpha2(1,t+1)=pavg(1,t+1);
    %end

    %RR1=[1; br(1,t); br(2,t)]*[1 , br(1,t), br(2,t)] + diag([0, wn(1,t)^2, wn(1,t)^2]);
    if ~all(isfinite(RR1))
        %warning = ('Matrix is Singular');
        break
    else
        if rank(RR1) < 3
            %warning = ('Matrix is Singular');
            break
        end
    end

    beta1(:,t+1)=beta1(:,t) ...
        +lambda*inv(RR1)*[1; p(1,t); p(2,t)]*(q1 -beta1(1,t)+beta1(2,t)*p(1,t)-beta1(3,t)*p(2,t));
    RR1 = RR1+lambda*( [ 1; p(1,t); p(2,t)]*[1, p(1,t), p(2,t)] - RR1);
    beta1(:,t+1)=max(0,beta1(:,t+1));

    %RR2=[1; br(2,t); br(1,t)]*[1 , br(2,t), br(1,t)] + diag([0, wn(2,t)^2, wn(2,t)^2]);
    if ~all(isfinite(RR2))
        %warning = ('Matrix is Singular');
        break
    else
        if rank(RR2) < 3
            %warning = ('Matrix is Singular');
            break
        end
    end

    beta2(:,t+1)=beta2(:,t) ...
        +lambda*inv(RR2)*[1; p(2,t); p(1,t)]*(q2 -beta2(1,t)+beta2(2,t)*p(2,t)-beta2(3,t)*p(1,t));
    RR2 = RR2+lambda*( [ 1; p(2,t); p(1,t)]*[1, p(2,t), p(1,t)] - RR2);
    beta2(:,t+1)=max(0,beta2(:,t+1));

end
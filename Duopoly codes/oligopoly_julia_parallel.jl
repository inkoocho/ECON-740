using LinearAlgebra, Plots

AA = 1;
BB = 1;
CC = 0.7;

# NE
p_nash = AA / (2*BB-CC);
p_collude=AA /(2*(BB-CC));

N = 30000;
lambda = 0.0025;
sigmawn = 0.0065;
include("oligopoly_julia_function.jl");
simul = 1;

result = zeros(Float64,simul,15,N);

@time Threads.@threads for i in 1:simul  
    result[i,:,:] = oligopoly_julia_function(N,lambda,sigmawn);
end

P = plot(layout = (3, 1), legendfontsize=8, size=(800, 1200), legend=:topleft)

for i in 1:simul
    p = result[i,:,:];
    plot!(P[1],p[1,1:end-1], alpha=1, label="Firm 1", linewidth=3)
    plot!(P[1],p[2,1:end-1], alpha=0.7, label="Firm 2", linewidth=2.5)
    plot!(P[1],p[3,1:end-1], alpha=0.4, label="Firm 3", linewidth=2)
    ylims!(P[1],p_nash-0.4, p_collude+0.2)
    plot!(P[1],[p_collude], seriestype="hline", label="Collusion Price")
    plot!(P[1],[p_nash], seriestype="hline", label="Nash Price")       
end


plot!(P[2],result[1,4,1:end-1], alpha=1, label="alpha0_12", linewidth=3)
plot!(P[2],result[1,6,1:end-1], alpha=0.7, label="alpha0_13", linewidth=2.5)

plot!(P[3],result[1,5,1:end-1], alpha=1, label="alpha1_12", linewidth=2)
plot!(P[3],result[1,7,1:end-1], alpha=0.7, label="alpha1_13", linewidth=1.5)
display(P)
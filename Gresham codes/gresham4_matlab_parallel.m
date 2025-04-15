pkg load parallel;
n=10; 
[pi,beta0,beta1,price] = pararrayfun(nproc, @gresham4_matlab_function,1:n, "Vectorized",true, "ChunksPerProc", 1, "UniformOutput", false); 
save gresham_result.mat pi beta0 beta1 price -mat7-binary
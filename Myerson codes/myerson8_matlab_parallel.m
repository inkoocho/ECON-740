pkg load statistics;
pkg load parallel;
tic
func_output = pararrayfun(nproc, @(n) myerson8_matlab_function(n),1:1000,"Vectorized",true, "ChunksPerProc", 10,"UniformOutput", false);
t4 = toc;
display(t4)
save myerson_result.mat func_output -mat7-binary
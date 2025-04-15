pkg load parallel;

sigmaf = linspace(0.0001,0.001,10);
sigmab = linspace(0.0001,0.001,10);
[sigmaf_grid, sigmab_grid] = meshgrid(sigmaf,sigmab);

tic
func_output = pararrayfun(100, @(sigmaf,sigmab) compete_fun(sigmaf,sigmab),sigmaf_grid,sigmab_grid);
t4 = toc;
display(t4)
save result_compete.mat func_output -mat7-binary
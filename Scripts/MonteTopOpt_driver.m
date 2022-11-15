%% Driver function to test the Layered Stochastic Topology optimization

% Generate mesh
nelemx = 100;
nelemy = 40;
box_mesh = box_mesh_generation_cantilever(nelemx,nelemy);

% run topology optimization
[xphy,E,E_C,Var_C] = abaqus_monte_1Dlayered_topCANTILIVER(box_mesh);
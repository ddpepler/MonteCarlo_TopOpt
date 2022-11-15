function [E] = abaqus_Rand_basis_layered(mesh,E_mean,E_std,corr_length)
% Generate a correlated random field using correlated guassian basis
% functions, and samples from a lognormal distribution

% number of layers
nl = mesh.nL;

% all the adhesion layer elements
ad_elem = mesh.elemnode((mesh.nelem+1):end,:);
nelem = size(ad_elem,1);
E = zeros(nelem,1);

%take the bottom line connections as the line elements for the 1D line for
% each layer
ad_line_elem = ad_elem(:,2:3);

for layer = 1:nl
    l_elem = find(mesh.elemlayerindex==layer);
    a = ad_line_elem(l_elem,:);
    [c,~,ic] = unique(a);
    meshline.elemnode = reshape(ic, size(a));
    meshline.coord = mesh.coord(c,:);
    [EL] = Random_field_basis_method(meshline,E_mean,E_std,corr_length);
    E(l_elem) = EL;
end
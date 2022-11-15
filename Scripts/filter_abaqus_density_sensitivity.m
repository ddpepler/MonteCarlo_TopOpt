function[dc_fac] = filter_abaqus_density_sensitivity(mesh,elem_area)
nelem = size(mesh.filter_fac,1);
weight_sum = sum(mesh.filter_fac,2);
dc_fac = mesh.filter_fac./repmat(weight_sum',nelem,1);
%dc_fac = diag(mesh.filter_fac)./(weight_sum);
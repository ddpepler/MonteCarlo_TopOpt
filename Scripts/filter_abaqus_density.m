function[xn] = filter_abaqus_density(mesh,x,elem_area)

weight_sum = mesh.filter_fac*elem_area;
xn = mesh.filter_fac*(x.*elem_area)./(weight_sum);
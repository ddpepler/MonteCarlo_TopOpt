function[mesh] = prepare_Sigmund_filter(mesh,rmin)
[nelem,~] = size(mesh.elemnode);
mp_temp = mesh.elem_midpoint(1,:);   
r_temp = sqrt((mesh.elem_midpoint(:,1)-mp_temp(1)).^2+(mesh.elem_midpoint(:,2)-mp_temp(2)).^2);    
num_points = length(find(rmin-r_temp'>0))*4;
fac_val = zeros(nelem,num_points);
fac_i = ones(nelem,num_points);
fac_j = ones(nelem,num_points);
for elem = 1:nelem
mp_temp = mesh.elem_midpoint(elem,:);   
r_temp = sqrt((mesh.elem_midpoint(:,1)-mp_temp(1)).^2+(mesh.elem_midpoint(:,2)-mp_temp(2)).^2);    
fac_temp = find(rmin-r_temp'>0);
fac_j(elem,1:length(fac_temp)) = fac_temp;
fac_i(elem,1:length(fac_temp)) = repmat(elem,[1,length(fac_temp)]);
fac_val(elem,1:length(fac_temp)) = rmin-r_temp(fac_temp)';
end
fac = sparse(fac_i,fac_j,fac_val);
mesh.filter_fac = fac;

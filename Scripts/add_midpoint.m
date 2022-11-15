function[mesh] = add_midpoint(mesh)

[nelem, ~] = size(mesh.elemnode);

midpoints = zeros(nelem,2);

for elem = 1:nelem
    quad = mesh.elemnode(elem,:);
    midpoints(elem,:) = [mean(mesh.coord(quad,1)),mean(mesh.coord(quad,2))];
end

mesh.elem_midpoint = midpoints;
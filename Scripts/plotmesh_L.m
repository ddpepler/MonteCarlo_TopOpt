function plotmesh_L(Mesh,density)

[nelem, nnode] = size(Mesh.elemnode);

meshx = zeros(4,2,nelem);
colormap(flipud(gray))
hold on
for elem = 1:nelem
   meshx(:,:,elem) = Mesh.coord(Mesh.elemnode(elem,1:4)',:); 
   fill(meshx(:,1,elem),meshx(:,2,elem),density(elem),'LineStyle','none')
% if mod(elem,200)==1
%    keyboard
% end
end
axis equal; axis tight
hold off
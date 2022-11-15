function plotmesh(Mesh,density)

[nelem, ~] = size(Mesh.elemnode);

meshx = zeros(4,2,nelem);
hold on
for elem = 1:nelem
   meshx(:,:,elem) = Mesh.coord(Mesh.elemnode(elem,1:4)',:); 
   x = 1-density(elem);
   if isnan(x)
   fill(meshx(:,1,elem),meshx(:,2,elem),[1,0,0],'LineStyle','none')    
   else
   %fill(meshx(:,1,elem),meshx(:,2,elem),[1,x,x],'LineStyle','none','FaceAlpha',0.5)
   fill(meshx(:,1,elem),meshx(:,2,elem),[x,x,x],'LineStyle','none')
   end

end
axis equal; axis tight; axis off
hold off
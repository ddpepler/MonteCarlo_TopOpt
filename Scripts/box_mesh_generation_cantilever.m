function[box_mesh] = box_mesh_generation_cantilever(nelemx,nelemy)
% nelemx= 20;
% nelemy = 50;
x = 0:1:nelemx;
y = 0:1:nelemy;
[X,Y] = meshgrid(x,y);
xy = [X(:), Y(:)];
nelem = nelemx*nelemy;
elems = zeros(nelem,4);
xc = 0:1:nelemx-1;
yc = nelemy:-1:1;
[XC,YC] = meshgrid(xc,yc);
for i = 1: nelem
   topleft = [XC(i), YC(i)];
   bottomleft = [XC(i), YC(i)-1];
   topright = [XC(i)+1, YC(i)];
   bottomright = [XC(i)+1, YC(i)-1];
   x1 = find(xy(:,1)==topleft(1));
   tlu = find(xy(x1,2) == topleft(2));
   tl = x1(tlu);
   x1 = find(xy(:,1)==topright(1));
   tru = find(xy(x1,2) == topright(2));
   tr = x1(tru);
   x1 = find(xy(:,1)==bottomleft(1));
   blu = find(xy(x1,2) == bottomleft(2));
   bl = x1(blu);
   x1 = find(xy(:,1)==bottomright(1));
   bru = find(xy(x1,2) == bottomright(2));
   br = x1(bru);
   elems(i,:) = [tl bl br tr]; 
end
box_mesh.coord = xy;
box_mesh.elemnode = elems;
% Tensile Specimen
% bottom_nodes= find(xy(:,2)==0);
% bc = bottom_nodes;
% Cantilever Beam
bc = find(xy(:,1)==0);
box_mesh.bc = bc;
%Tensile Specimen
Right_nodes= find(xy(:,1)== nelemx);
Top_of_right_nodes = find(xy(Right_nodes,2) == nelemy);
displaced_nodes = Right_nodes(Top_of_right_nodes);
% displaced_nodes = find(xy(:,1) == nelemx && xy(:,2) == 0);
box_mesh.ld = displaced_nodes;
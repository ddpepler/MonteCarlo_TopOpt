function[mesh] = add_adhesionlayers_quad(mesh,lt)
%version #1 - adds an adhesion layer elements for a given mesh. Currently
%set up for horizontal layers
%version #2 - created a mapping from original mesh to new mesh for density
%distribution
%version #3 - added a layer index and layer count to the mesh structure
%version #5 - changed to quadtratic elements

[nelem, nshp] = size(mesh.elemnode);
[nnode] = size(mesh.coord,1);
new_nelem = nelem*2;
new_nodes = new_nelem*5;
new_coord = zeros(nnode+new_nodes,2);
new_elemnode = zeros(nelem+new_nelem,nshp);
new_coord(1:nnode,:) = mesh.coord;
topgroup = zeros(nelem+new_nelem,1);
%layer_thickness = 0.01; %percentage of regular element
layer_thickness = lt;
count = 1;

elem1node = mesh.elemnode(1,:);
%% mesh spacing

% X
%h = abs(mesh.coord(elem1node(2),1)-mesh.coord(elem1node(3),1));

% Y
h = abs(mesh.coord(elem1node(1),2)-mesh.coord(elem1node(2),2));

eps = h*layer_thickness;


%% add nodes and elements

for elem = 1:nelem
    quad = mesh.elemnode(elem,:);
    bottom_nodes = [mesh.coord(quad(2:3),:);mesh.coord(quad(6),:)];
    if bottom_nodes(1,2) == 0
    else
        if mesh.coord(quad(2),1) == 0
           new_coord(nnode+count,:) = [bottom_nodes(1,1), bottom_nodes(1,2)+eps];
           count = count+1;
           new_coord(nnode+count,:) = [bottom_nodes(2,1), bottom_nodes(2,2)+eps];
           count = count+1;
           new_coord(nnode+count,:) = [bottom_nodes(3,1), bottom_nodes(3,2)+eps];
           count = count+1;
           new_coord(nnode+count,:) = [bottom_nodes(1,1), bottom_nodes(1,2)+eps/2];
           count = count+1;
           new_coord(nnode+count,:) = [bottom_nodes(2,1), bottom_nodes(2,2)+eps/2];
           count = count+1;
           new_coord(quad(5),2) = mesh.coord(quad(1),2)-(mesh.coord(quad(1),2)-mesh.coord(quad(2),2)-eps)/2;
           new_coord(quad(7),2) = mesh.coord(quad(4),2)-(mesh.coord(quad(4),2)-mesh.coord(quad(3),2)-eps)/2;
      
        else
            
          new_coord(nnode+count,:) = [bottom_nodes(2,1), bottom_nodes(2,2)+eps];  
          count = count+1; 
          new_coord(nnode+count,:) = [bottom_nodes(3,1), bottom_nodes(3,2)+eps];
          count = count+1;
          new_coord(nnode+count,:) = [bottom_nodes(2,1), bottom_nodes(2,2)+eps/2];
          count = count+1;
          new_coord(quad(5),2) = mesh.coord(quad(1),2)-(mesh.coord(quad(1),2)-mesh.coord(quad(2),2)-eps)/2;
          new_coord(quad(7),2) = mesh.coord(quad(4),2)-(mesh.coord(quad(4),2)-mesh.coord(quad(3),2)-eps)/2;

        end
    end
end

new_coord = new_coord(1:nnode+count-1,:);

%% create new elements
count = 1;
for elem = 1:nelem
    quad = mesh.elemnode(elem,:);
    bottom_nodes = [mesh.coord(quad(2:3),:);mesh.coord(quad(6),:)];
    if bottom_nodes(1,2) == 0
        new_elemnode(elem,:) = quad;
        topgroup(elem) = elem;
    else
        nbn1 = find(new_coord(:,2)==bottom_nodes(1,2)+eps);
        bn1 = find(new_coord(nbn1,1) == bottom_nodes(1,1));
        nbn2 = find(new_coord(:,2)==bottom_nodes(2,2)+eps);
        bn2 = find(new_coord(nbn2,1) == bottom_nodes(2,1));
        nbn3 = find(new_coord(:,2)==bottom_nodes(3,2)+eps);
        bn3 = find(new_coord(nbn3,1) == bottom_nodes(3,1));
        nbn4 = find(new_coord(:,2)==bottom_nodes(1,2)+eps/2);
        bn4 = find(new_coord(nbn4,1) == bottom_nodes(1,1));
        nbn5 = find(new_coord(:,2)==bottom_nodes(2,2)+eps/2);
        bn5 = find(new_coord(nbn5,1) == bottom_nodes(2,1));
        bn1 = nbn1(bn1(1));
        bn2 = nbn2(bn2(1));
        bn3 = nbn3(bn3(1));
        bn4 = nbn4(bn4(1));
        bn5 = nbn5(bn5(1));
        new_elemnode(elem,:) = [quad(1), bn1, bn2, quad(4), quad(5), bn3, quad(7), quad(8)];
        new_elemnode(nelem+count,:) = [bn1,quad(2), quad(3),bn2, bn4, quad(6), bn5, bn3];
        topgroup(elem) = elem;
        topgroup(nelem+count) = elem;            
        count = count+1;
    end

end

new_elemnode = new_elemnode(1:nelem+count-1,:);
topgroup = topgroup(1:nelem+count-1,1);
mesh.topgroup = topgroup;
mesh.coord = new_coord;
mesh.elemnode = new_elemnode;

% Version #3 add adhesion layer numbering

nal = count-1;
elem_al = new_elemnode(nelem+1:nelem+nal,:);
midpoint_al = zeros(nal,2);
for i = 1:nal
   quad = elem_al(i,:);
   midpoint_al(i,:) = [mean(mesh.coord(quad,1)),mean(mesh.coord(quad,2))]; 
end

[C, ~ ,ia] = unique(midpoint_al(:,2));

nL = size(C,1);
mesh.elemlayerindex = ia;
mesh.nL = nL;
mesh.nelem = nelem;



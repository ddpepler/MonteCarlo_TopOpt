function[mesh] = add_adhesionlayers_4(mesh,lt)
%version #1 - adds an adhesion layer elements for a given mesh. Currently
%set up for horizontal layers
%version #2 - created a mapping from original mesh to new mesh for density
%distribution
%version #3 - added a layer index and layer count to the mesh structure

[nelem, nshp] = size(mesh.elemnode);
[nnode] = size(mesh.coord,1);
new_nelem = nelem*2-1;
new_nodes = new_nelem+1;
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


%% add nodes

for elem = 1:nelem
    quad = mesh.elemnode(elem,:);
    bottom_nodes = mesh.coord(quad(2:3),:);
    
    if bottom_nodes(1,2) == 0
    else
        if mesh.coord(quad(2),1) == 0
           new_coord(nnode+count,:) = [bottom_nodes(1,1), bottom_nodes(1,2)+eps];
           count = count+1;
           new_coord(nnode+count,:) = [bottom_nodes(2,1), bottom_nodes(2,2)+eps];
           count = count+1;
      
        else
            
          new_coord(nnode+count,:) = [bottom_nodes(2,1), bottom_nodes(2,2)+eps];  
          count = count+1;  

        end
    end
end

new_coord = new_coord(1:nnode+count-1,:);

%% create new elements
count = 1;
for elem = 1:nelem
    quad = mesh.elemnode(elem,:);
    bottom_nodes = mesh.coord(quad(2:3),:);
    if bottom_nodes(1,2) == 0
        new_elemnode(elem,:) = quad;
	topgroup(elem) = elem;
    else
        nbn1 = find(new_coord(:,2)==bottom_nodes(1,2)+eps);
        bn1 = find(new_coord(nbn1,1) == bottom_nodes(1,1));
        nbn2 = find(new_coord(:,2)==bottom_nodes(2,2)+eps);
        bn2 = find(new_coord(nbn2,1) == bottom_nodes(2,1));
        bn1 = nbn1(bn1(1));
        bn2 = nbn2(bn2(1));
        new_elemnode(elem,:) = [quad(1), bn1, bn2, quad(4)];
        new_elemnode(nelem+count,:) = [bn1,quad(2), quad(3),bn2];
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
midpoint_bul = zeros(nelem,2);
for j = 1:nelem
    quad = mesh.elemnode(j,:);
    midpoint_bul(j,:) = [mean(mesh.coord(quad,1)),mean(mesh.coord(quad,2))];
end

[C, ~ ,ia] = unique(midpoint_al(:,2));
[~,~,ib] = unique(midpoint_bul(:,2));

nL = size(C,1);
mesh.elemlayerindex = ia;
mesh.bul_elemlayerindex = ib;
mesh.nL = nL;
mesh.nelem = nelem;



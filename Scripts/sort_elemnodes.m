function [mesh] = sort_elemnodes(mesh)

[nelem, nshp] = size(mesh.elemnode);
if nshp == 4
    for elem = 1:nelem
       quad = mesh.elemnode(elem,:); 
       maxx = max(mesh.coord(quad,1));
       maxy = max(mesh.coord(quad,2));
       minx = min(mesh.coord(quad,1));
       miny = min(mesh.coord(quad,2));

       nnq1 = find(mesh.coord(quad,1)==minx);
       nq1 = find(mesh.coord(quad(nnq1),2)==maxy);
       nnq2 = find(mesh.coord(quad,1)==minx);
       nq2 = find(mesh.coord(quad(nnq2),2)==miny);
       nnq3 = find(mesh.coord(quad,1)==maxx);
       nq3 = find(mesh.coord(quad(nnq3),2)==miny);
       nnq4 = find(mesh.coord(quad,1)==maxx);
       nq4 = find(mesh.coord(quad(nnq4),2)==maxy);
       mesh.elemnode(elem,:) = [quad(nnq1(nq1)),quad(nnq2(nq2)),quad(nnq3(nq3)),quad(nnq4(nq4))];
    end
elseif nshp == 8
    for elem = 1:nelem
       quad = mesh.elemnode(elem,:);
       maxx = max(mesh.coord(quad,1));
       maxy = max(mesh.coord(quad,2));
       minx = min(mesh.coord(quad,1));
       miny = min(mesh.coord(quad,2));
       [~, ix] = sort(mesh.coord(quad,1));
       [~, iy] = sort(mesh.coord(quad,2));
       ln = quad(ix(1:3));
       mn1 = quad(ix(4:5));
       rn = quad(ix(6:end));
       md2 = quad(iy(4:5));
       q1 = ln(find(mesh.coord(ln,2) == maxy));
       q2 = ln(find(mesh.coord(ln,2) == miny));
       q3 = rn(find(mesh.coord(rn,2) == miny));
       q4 = rn(find(mesh.coord(rn,2) == maxy));
       q5 = md2(find(mesh.coord(md2,1) == minx));
       q6 = mn1(find(mesh.coord(mn1,2) == miny));
       q7 = md2(find(mesh.coord(md2,1) == maxx));
       q8 = mn1(find(mesh.coord(mn1,2) == maxy));
       mesh.elemnode(elem,:) = [q1, q2 , q3, q4, q5, q6, q7, q8];
    end
end

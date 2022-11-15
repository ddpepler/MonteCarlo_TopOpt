function [shp,shpx] = monomial_sq(x,elem_order)

x1 = x(:,1);
x2 = x(:,2);
np = size(x,1);
oo = ones(np,1);
zz = zeros(np,1);
switch elem_order
    case 1
shp = [x2, oo, x1, x1.*x2 ];
shpx(:,:,1) = [zz, zz, oo, x2];
shpx(:,:,2) = [oo, zz, zz, x1];
    case 2
shp = [x2, oo, x1, x1.*x2, x2.^2, x1.^2, x1.*x2.^2,x2.*x1.^2 ];
shpx(:,:,1) = [zz, zz, oo, x2, zz, 2*x1, x2.^2, 2*x1.*x2];
shpx(:,:,2) = [oo, zz, zz, x1, 2*x2, zz, 2*x1.*x2, x1.^2];
    case 3
shp = [x2, oo, x1, x1.*x2, x2.^2, x1.^2, x1.*x2.^2,x2.*x1.^2, x2.^3, x1.^3, x1.*x2.^3, x2.*x1.^3];
shpx(:,:,1) = [zz, zz, oo, x2, zz, 2*x1, x2.^2, 2*x1.*x2, zz, 3*x1.^2, x2.^3, 3*x2.*x1.^2];
shpx(:,:,2) = [oo, zz, zz, x1, 2*x2, zz, 2*x1.*x2, x1.^2, 3*x2.^2, zz, 3*x1.*x2.^2, x1.^3];        

end
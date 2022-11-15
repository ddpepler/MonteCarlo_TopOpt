% Shape function on reference element 2D quad

function ref = make_ref_square(elem_order)
if nargin < 1
    elem_order = 1;
end
% quadrature set up
detJ = 1/4;
a = 1/sqrt(3);
ref.xq = [(1-a)/2 (1-a)/2; (1-a)/2 (1+a)/2; (1+a)/2 (1-a)/2; (1+a)/2 (1+a)/2];
ref.wq = detJ*[1; 1; 1; 1];
switch elem_order
    case 1
% interpl nodes
xnodes = [0 1; 0 0; 1 0; 1 1];
    case 2
xnodes = [0 1; 0 0; 1 0; 1 1; 0 0.5; 0.5 0; 1 0.5; 0.5 1]; 
    case 3
xnodes = [0 1; 0 0; 1 0; 1 1; 0 1/3; 1/3 0; 1 1/3; 1/3 1; 0 2/3; 2/3 0; 1 2/3; 2/3 1];        
end
inv_coeff = monomial_sq(xnodes,elem_order);
[psi, psix] = monomial_sq(ref.xq,elem_order);

ref.shp = psi/inv_coeff;
ref.shpx(:,:,1) = psix(:,:,1)/inv_coeff;
ref.shpx(:,:,2) = psix(:,:,2)/inv_coeff;
end
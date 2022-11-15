function[amat] = stoch_sq_fe(mesh,Kin)
%Inputs:
%mesh - The mesh structure from abaqus (element node number relations, node
%coordinates, load nodes, boundary condition nodes)
%Kin - The nth term in the Ki expansion for the Young's modulus (nelem X 1
%vector)

%material parameters
nu = 0.3;
lam = nu/(1-nu^2);
mu = 1/(2*(1+nu));


[nelem, nshp] = size(mesh.elemnode);
if nshp == 8
    p = 2;
else
    p = 1;
end
ref = make_ref_square(p);
ldof{1} = 1:2:2*nshp;
ldof{2} = 2:2:2*nshp;
nldof = 2*nshp;
amat = zeros(nldof,nldof,nelem);

for elem = 1:nelem
    
    quad = mesh.elemnode(elem,:)';
    
    %jacobian calculation
    xl = mesh.coord(quad,:);
    jacq11 = ref.shpx(:,:,1)*xl(:,1);
    jacq12 = ref.shpx(:,:,2)*xl(:,1);
    jacq21 = ref.shpx(:,:,1)*xl(:,2);
    jacq22 = ref.shpx(:,:,2)*xl(:,2);
    detjq = jacq11.*jacq22 - jacq12.*jacq21;
    wqJ = ref.wq.*detjq;
    %inverse jacobian
    itjac11 = jacq22./detjq;
    itjac12 = -jacq21./detjq;
    itjac21 = -jacq12./detjq;
    itjac22 = jacq11./detjq;
    phixq11 = bsxfun(@times,ref.shpx(:,:,1), repmat(itjac11,1,nshp));
    phixq21 = bsxfun(@times,ref.shpx(:,:,1), repmat(itjac21,1,nshp));
    phixq12 = bsxfun(@times,ref.shpx(:,:,2), repmat(itjac12,1,nshp));
    phixq22 = bsxfun(@times,ref.shpx(:,:,2), repmat(itjac22,1,nshp));
    phixq1 = phixq11+phixq12;
    phixq2 = phixq21+phixq22;

    %compute local sitffness matrices
    aloc = zeros(nldof,nldof);
    aloc(ldof{1},ldof{1}) = mu*phixq2'*diag(wqJ)*phixq2 + 2*mu*phixq1'*diag(wqJ)*phixq1 + lam*phixq1'*diag(wqJ)*phixq1; 
    aloc(ldof{1},ldof{2}) = mu*phixq2'*diag(wqJ)*phixq1 + lam*phixq1'*diag(wqJ)*phixq2;
    aloc(ldof{2},ldof{1}) = mu*phixq1'*diag(wqJ)*phixq2 + lam*phixq2'*diag(wqJ)*phixq1;
    aloc(ldof{2},ldof{2}) = mu*phixq1'*diag(wqJ)*phixq1 + 2*mu*phixq2'*diag(wqJ)*phixq2 + lam*phixq2'*diag(wqJ)*phixq2;  
    amat(:,:,elem) = Kin(elem)*aloc;
end
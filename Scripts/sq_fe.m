function[amat,imat,jmat,elem_area] = sq_fe(mesh)

%material parameters
nu = 0.3;
lam = nu/(1-nu^2);
mu = 1/(2*(1+nu));
E = 71000;


[nelem, nshp] = size(mesh.elemnode);
if nshp == 12
    p =3;
elseif nshp == 8
    p = 2;
else
    p = 1;
end
ref = make_ref_square(p);
ldof{1} = 1:2:2*nshp;
ldof{2} = 2:2:2*nshp;
nldof = 2*nshp;
amat = zeros(nldof,nldof,nelem);
imat = zeros(nldof,nldof,nelem);
jmat = zeros(nldof,nldof,nelem);
elem_area = zeros(nelem,1);
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
    elem_area(elem) = sum(wqJ);
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
    iloc = zeros(nldof,nldof);
    jloc = zeros(nldof,nldof);
    aloc(ldof{1},ldof{1}) = mu*phixq2'*diag(wqJ)*phixq2 + 2*mu*phixq1'*diag(wqJ)*phixq1 + lam*phixq1'*diag(wqJ)*phixq1; 
    aloc(ldof{1},ldof{2}) = mu*phixq2'*diag(wqJ)*phixq1 + lam*phixq1'*diag(wqJ)*phixq2;
    aloc(ldof{2},ldof{1}) = mu*phixq1'*diag(wqJ)*phixq2 + lam*phixq2'*diag(wqJ)*phixq1;
    aloc(ldof{2},ldof{2}) = mu*phixq1'*diag(wqJ)*phixq1 + 2*mu*phixq2'*diag(wqJ)*phixq2 + lam*phixq2'*diag(wqJ)*phixq2;
    iloc(ldof{1},ldof{1}) = repmat(2*quad-1,[1,nshp]);
    iloc(ldof{1},ldof{2}) = repmat(2*quad-1,[1,nshp]);
    iloc(ldof{2},ldof{1}) = repmat(2*quad,[1,nshp]);
    iloc(ldof{2},ldof{2}) = repmat(2*quad,[1,nshp]);
    jloc(ldof{1},ldof{1}) = repmat(2*quad'-1,[nshp, 1]);
    jloc(ldof{1},ldof{2}) = repmat(2*quad',[nshp, 1]);
    jloc(ldof{2},ldof{1}) = repmat(2*quad'-1,[nshp, 1]);
    jloc(ldof{2},ldof{2}) = repmat(2*quad',[nshp, 1]);
        
    amat(:,:,elem) = E*aloc;
    imat(:,:,elem) = iloc;
    jmat(:,:,elem) = jloc;
end

function[xphy,E,E_C,Var_C] = abaqus_monte_1Dlayered_topCANTILIVER(mesh)

% performs topology optimization on a mesh structure from abaqus
[nelem_map,nshp] = size(mesh.elemnode);
%check element node numering and fix any that are numbered incoorectly
[mesh] = sort_elemnodes(mesh);
%add adhesion layers
mesh = add_midpoint(mesh);
lt = 0.2;
rmin = 1.5*2;
mesh = prepare_Sigmund_filter(mesh,rmin);
%add adhesion layers
%linear elements
if nshp == 4
[mesh2] = add_adhesionlayers_4(mesh,lt);

%quadratic elements
elseif nshp == 8
[mesh2] = add_adhesionlayers_quad(mesh,lt);
else
    error('unsupported mesh type')
end
mesh2 = add_midpoint(mesh2);
%topology optimization parameters
volfrac = 1;
penal = 3;
kappa = 1;
[nelem,nshp] = size(mesh2.elemnode);
nnode = size(mesh2.coord,1);
x_map = mesh2.topgroup;
nf = 500; %number of random fields considered in the Monte-Carlo
%% 1D Ransom Fields from Experiments
[E] = abaqus_layered_EfieldGen(mesh2,nf);
E = E*71000;
%% 1D Random Variables
% % Generate random fields
% E = ones(nelem,nf);
% for i = 1:nf
% El = abaqus_Rand_basis_layered(mesh2,E_mean,E_std,corr_length);
% E((nelem_map+1):end,i) = El;
% end
% E = E*71000;
%% uniform decrease in layer strength
%E((nelem_map+1):end,1) = 0.1;
%% graded layer weakening
% x1 = 0.2;
% x2 = 0.5;
% El = linspace(x1,x2,mesh2.nL);
% for i = 1:mesh2.nL
% 	I = find(mesh2.elemlayerindex==i);
% 	I = I + (nelem_map); 
% 	E(I) = El(i);
% end
%% Random E value for each layer
% El = 0.1*rand(mesh2.nL,nf)+0.5;
% for i = 1:mesh2.nL
% 	I = find(mesh2.elemlayerindex==i);
% 	I = I + (nelem_map); 
%     ni = size(I,1);
% 	E(I,:) = repmat(El(i,:),ni,1);
% end
xL = ones(nelem,1)*volfrac;
x = ones(nelem_map,1)*volfrac;
xold1 = x;
xold2 = x;
%create local stiffness matrices for FEA
[amat,imat,~,elem_area] = sq_fe(mesh);

[dc_fac] = filter_abaqus_density_sensitivity(mesh,elem_area);
%assemble Ki matricies
smatL = zeros(2*nshp,2*nshp,nelem,nf);
for i = 1:nf
[smatL(:,:,:,i)] = stoch_sq_fe(mesh2,E(:,i));
end
[~,imatL,jmatL,~] = sq_fe(mesh2);
total_volume = sum(elem_area);
F = zeros(2*nnode,1);
%U = zeros(2*nnode,1);
F(2*mesh.ld) = -200/length(mesh.ld); %Force in newtons, up is positive
%fixeddof = [2*mesh.bc1-1; 2*mesh.bc2];
fixeddof = [2*mesh.bc-1; 2*mesh.bc ];
alldof = 1:2*nnode;
freedofs = setdiff(alldof,fixeddof);

c_max = 200;
change = 1;
loop = 0;
m_mma =1; %number of constraints
low = 0;
upp = 0;
beta = 0.01;
[xphy] = filter_abaqus_density(mesh,x,elem_area);
[dc_fac_H,xphy] = filter_abaqus_heaviside(beta,xphy);
xphyL = xphy(x_map);
while change > 0.001

    loop = loop +1;
    C_i = zeros(nf,1);
    dc_i = zeros(nelem_map,nf);
    %% FEA
    %amat_it = zeros(size(amatL));
    smat_it = zeros(size(smatL));
    %dc = zeros(nelem_map,1);
    for i = 1: nelem
        smat_it(:,:,i,:) = smatL(:,:,i,:)*xphyL(i)^penal;
    end
    parfor (j =1:nf) %parallel loop used, can switch to a regular for loop if uses to much CPU or RAM
        U = zeros(2*nnode,1);
        amat_t = smat_it(:,:,:,j);
        K = sparse(imatL(:),jmatL(:),amat_t(:));
        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
        [m,n] = size(K);
        C_i(j) = U'*K*U;
        %% sensitivities
        for i = 1: nelem_map
           U_temp = zeros(2*nnode,1);
           loc_index = find(x_map==i);
           quad = imatL(:,1,loc_index);
           quad = unique(quad);
           U_temp(quad) = U(quad);
           imat_loc = imatL(:,:,loc_index);
           jmat_loc = jmatL(:,:,loc_index);
           amat_loc = smatL(:,:,loc_index,j);
           K_loc = sparse(imat_loc(:),jmat_loc(:),amat_loc(:),m,n);
           dc_i(i,j) = -penal*xphy(i)^(penal-1)*U_temp'*K_loc*U_temp;
        end
    end
    
    E_C = (1/nf)*sum(C_i);
    Var_C = (1/nf)*sum((C_i-E_C).^2);
    if Var_C < 1e-12 % condition for when variance is 0
        dc = zeros(nelem_map,1);
        for m = 1:nf 
        dc = dc + dc_i(:,m)./nf;  %limit for sensitivity as Var_C approaches 0
        end    
    else
        dc = zeros(nelem_map,1);
        for m = 1:nf 
        dc = dc + (1/nf+kappa*((C_i(m)-E_C))/(sqrt(Var_C)*nf))*dc_i(:,m); %analytic expression for sensitivity 
        end
    end
    dc = dc;
     C_hat = E_C + kappa*sqrt(Var_C);
    %% filter function
    %if loop <20.5
    %[dc] = filter_abaqus(mesh,dc);
%     [dc] = filter_abaqus_2(mesh,dc);
    
    %end
    dc = (dc_fac.*dc_fac_H)*dc;
    %% Optimality function
    %     %min compliance
%     fval = sum(x.*elem_area)-total_volume*volfrac;
%     dfdx = elem_area';
%     [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m_mma,nelem_map,loop,x ... 
%       ,0.001,1,xold1,xold2,c,dc,fval,dfdx,low,upp,1,0,1e6,1,1,0.5);
%   xold2 = xold1;
%   xold1 = x;
%   x = xmma;
%   xL = x(x_map);
    fval = sum(x.*elem_area);
    dfdx = elem_area;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m_mma,nelem_map,loop,x ... 
      ,0.001,1,xold1,xold2,fval,dfdx,C_hat-c_max,dc',low,upp,1,0,1e6,1,0.125,0.125);
  xold2 = xold1;
  xold1 = x;
  x = xmma;
  
  [xphy] = filter_abaqus_density(mesh,x,elem_area);
  [dc_fac_H,xphy] = filter_abaqus_heaviside(beta,xphy);
  if loop < 100
      if mod(loop,20)==0
        beta = beta*2;
      end
  end
  xphyL = xphy(x_map);
  
  change = min(sum(abs(x-xold1).*elem_area)/(total_volume*volfrac),sum(abs(x-xold2).*elem_area)/(total_volume*volfrac));

    
%   change1 = min(max(max(abs(x-xold1))),max(max(abs(x-xold2))));
%   total_change = [total_change c];
%   mean_new = mean(total_change);
%   %change = min(abs(mean_new-mean_old)/mean_new,change1);
%   change = change1;
%   mean_old = mean_new;
  
  
  %% print results
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',C_hat) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x.*elem_area))/total_volume) ...
        ' ch.: ' sprintf('%6.3f',change )])
    
  %% plot iteration
  %plotmesh(mesh,x)
  if loop > 1000
      print('max iterations reached')
      break
  end
  
end
  %% plot final solution
  figure(1)
  plotmesh(mesh,x)
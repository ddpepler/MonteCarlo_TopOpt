function[E] = abaqus_layered_EfieldGen(mesh,nf)

%% random field characteristics for fusion material
RFinput.Type = 'Lognormal';
% experimental
% TheCOV = 0.09 ; %covariance in terms of mean
% TheMean = 0.89; %mean of random field
% modifid
TheCOV = 0.2 ; %covariance in terms of mean
TheMean = 0.4; %mean of random field
if isequal(RFinput.Type ,'Lognormal') 
   RFinput.LNMean = TheMean  ;
   RFinput.LNStdv = TheCOV * RFinput.LNMean ; 
end
if isequal(RFinput.Type ,'Gaussian') 
   RFinput.Mean = TheMean  ;
   RFinput.Stdv = TheCOV * RFinput.Mean ;
end
RFinput.CorrType = 'exp' ; 
RFinput.CorrLength = 1000; %correlation length
RFinput.DiscScheme = 'KL'; %random field discretization scheme
RFinput.OrderExp = 4; %number of terms taken in eigen decomposition
max_coord = max(mesh.coord(:,1));
min_coord = min(mesh.coord(:,1));
RFinput.Domain = {min_coord, max_coord};
RF_ad = DiscRandomField([],RFinput);


%% random field characteristics for bulk material
% TheCOV = 0.068 ; %covariance in terms of mean
% TheMean = 1; %mean of random field
TheCOV = 0.2 ; %covariance in terms of mean
TheMean = 0.5; %mean of random field
if isequal(RFinput.Type ,'Lognormal') 
   RFinput.LNMean = TheMean  ;
   RFinput.LNStdv = TheCOV * RFinput.LNMean ; 
end
if isequal(RFinput.Type ,'Gaussian') 
   RFinput.Mean = TheMean  ;
   RFinput.Stdv = TheCOV * RFinput.Mean ;
end
% RFinput.CorrLength = 6.25;
RF_bul = DiscRandomField([],RFinput);

%% building the random field
nl = mesh.nL;
nelem = mesh.nelem;
E = zeros(size(mesh.elemnode,1),nf);
% deterministic bulk layer
% E(1:mesh.nelem,:) = 1;
for flds = 1:nf
%     random field 
    for bul_layer = 1:nelem
        l_elem = find(mesh.bul_elemlayerindex==bul_layer);
        xl = mesh.elem_midpoint(l_elem,1);
        basisl = zeros(length(xl),1);
        Xi = randn([RFinput.OrderExp,1]);
        for i = 1:length(xl)
            basisl(i) = EvalRandomField(RF_bul,xl(i),Xi);
        end
        E(l_elem,flds) = basisl;
    end
% deterministic bulk layer
    for layer = 1:nl
        l_elem = find(mesh.elemlayerindex==layer);
        xl = mesh.elem_midpoint(mesh.nelem+l_elem,1);
        basisl = zeros(length(xl),1);
        Xi = randn([RFinput.OrderExp,1]);
        for i = 1:length(xl)
%             Xi = randn([RFinput.OrderExp,1]);
            basisl(i) = EvalRandomField(RF_ad,xl(i),Xi);
        end
        E(mesh.nelem+l_elem,flds) = basisl;
    end
end

function RF = DiscGaussianRandomField(COORD,RFinput)
%  
%      RF = DiscGaussianRandomField(COORD, RFinput)
%      Compute data for gaussian random field discretization
%  
%  

%%------------------------------------------------------------%%
%            Common initialization for all discretization methods
%%------------------------------------------------------------%%

RF.Type = RFinput.Type ;
RF.Mean = RFinput.Mean ;      % homogeneous random fields are considered
RF.Stdv = RFinput.Stdv ;      

if isfield(RFinput, 'Domain')
   RF.Domain = RFinput.Domain;
   xmin = RF.Domain{1};
   xmax = RF.Domain{2};
else  %    computing automatically the size of discretization domain
   xmin = min(COORD);            
   xmax = max(COORD);
   RF.Domain = { xmin , xmax};
end

%    assigning input data
RF.CorrType   = RFinput.CorrType;
RF.CorrLength = RFinput.CorrLength;
RF.DiscScheme = RFinput.DiscScheme;
RF.OrderExp   = RFinput.OrderExp;
RF.Outcome    = zeros([RF.OrderExp , 1]);

OrderGaussExp = RF.OrderExp;

switch RF.DiscScheme
case 'KL'
   %%------------------------------------------------------------%%
   %            Karhunen Loeve expansion
   %%------------------------------------------------------------%%
   switch RF.CorrType
   case 'exp'
      %  SizeDomain is [ax , ay] for 2D problem such that
      %  the KL eigenvalue problem is solved over [-ax , ax] x [-ay , ay]
      SizeDomain = (xmax-xmin)/2;
      RF.Translation =(xmax + xmin)/2;
      %    computing basis translation vector T
      %    if the structure has arbitrary coordinates, the KL problem
      %    is solve after translation
      %    X_ref = X_user - T
      %    In computation phi(X_user) is obtained by EvalBasis(X_ref)
      %
      
      
      %    computing data for eigenvalues and eigenfunctions
      RFtmp = KLExpansion(SizeDomain , RF.CorrLength , OrderGaussExp);
      RF.EVPx = RFtmp.EVPx ;
      RF.Eigs = RFtmp.Eigs ;
      if (length(RF.CorrLength) ==2 ) 
         RF.EVPx = RFtmp.EVPx ;
         RF.EVPy = RFtmp.EVPy ;
         RF.Eigs = RFtmp.Eigs ;
         RF.EVP2pos = RFtmp.EVP2pos ;
      end;
      
   otherwise
      fprintf('\n\n   NB : The Karhunen Loeve Expansion is available\n');
      fprintf('          only for exponential-type correlation function\n\n');
      fprintf('   ------>  Program goes on with RF.CorrType = ''exp''\n');
      RFinput.CorrType = 'exp';
      RF = DiscGaussianRandomField(COORD,RFinput);
      
   end
   
case 'EOLE'
   %%------------------------------------------------------------%%
   %            EOLE expansion of 2D fields
   %%------------------------------------------------------------%%
   if isfield(RFinput, 'Npts')
      RF.Npts = RFinput.Npts ;	
   else
      fprintf('\n\n NB : Grid not specified for EOLE expansion \n');
      fprintf('      Uniform grid chosen with 5 points in the correlation length\n\n');
      
      switch length(RF.CorrLength)
      case 1
         RF.Npts = 1 + 5 * ceil((xmax(1) -xmin(1))/RF.CorrLength);
         fprintf('      RF.Npts = %3.0d\n\n',RF.Npts);
      case 2
         RF.Npts(1) = 1 + 5 * ceil((xmax(1) -xmin(1))/RF.CorrLength(1));
         RF.Npts(2) = 1 + 5 * ceil((xmax(2) -xmin(2))/RF.CorrLength(2));
         fprintf('      RF.Npts = [%3.0d , %3.0d]\n\n',RF.Npts(1),RF.Npts(2));
      end
   end
   
   switch length(RF.CorrLength)
   case 1
      VV = (xmin(1) : (xmax(1) -xmin(1))/(RF.Npts(1)-1) : xmax(1))';
   case 2
      [A , B ] = meshgrid(xmin(1) : (xmax(1) -xmin(1))/(RF.Npts(1)-1) : xmax(1) , ...
         xmin(2) : (xmax(2) -xmin(2))/(RF.Npts(2)-1) : xmax(2));
      
      Asize = size(A);
      NbRFnodes = Asize(1) * Asize(2);
      VV = cat(2,reshape(A',NbRFnodes,1),reshape(B',NbRFnodes,1));
   otherwise
   end;
   
   CorrMat = zeros(length(VV));
   for i = 1 : length(VV)
      for j = i : length(VV)
         CorrMat(i,j) = CorrFunEval(RF.CorrType , VV(i,:) , VV(j,:) , ...
            RF.CorrLength);
         CorrMat(j,i) = CorrMat(i,j);
      end;
   end;
   %CovMat;
   options.disp = 0;
   %  Verifies CovMat*Phi =  Phi*Theta
   [Phi , Theta] = eigs(CorrMat,eye(size(CorrMat)),OrderGaussExp);
   RF.COORD = VV;
   RF.Phi = Phi;              % eigenvectors of correlation/and
   % covariance matrix 
   RF.Eigs = diag(Theta);  % eigenvalues of CorrMatrix order from
   % the largest 
   
   
case 'OSE'
   SizeDomain = (xmax-xmin)/2;                 % a 
   RF.Translation =(xmax + xmin)/2;            % T
   switch length(RF.CorrLength)
   case 1
      %	if (OrderGaussExp <= 10)
      %	  NPG = max(OrderGaussExp,1) ;
      %	else
      %	  NPG = 16;
      %	end
      
      %   Maximal accuracy in integral computation 
      NPG  = 16 ;
      [xIP , WGH ] = GaussPoints(NPG);
      LP =  LegendrePolynomials(OrderGaussExp,xIP);
      
      %        Compute the correlation between all Gauss Points
      C = zeros(NPG);
      for i = 1 : NPG
         for j = i : NPG
            C(i,j) = CorrFunEval(RF.CorrType , SizeDomain * xIP(i), ...
               SizeDomain * xIP(j), ...
               RF.CorrLength);
            C(j,i) = C(i,j);
         end;
      end;
      %        Compute CovChi(k,l) = correlation matrix of Chi-variables
      %        !!!! to get the covariance  amtrix, result should be
      %        multiplied by Stdv^2.
      CovChi = zeros(OrderGaussExp);
      for k = 1 : OrderGaussExp
         for l = k : OrderGaussExp
            for i = 1 : NPG
               for j = 1 : NPG
                  CovChi(k,l) = CovChi(k,l) + ...
                     WGH(i) * WGH(j) * C(i,j) *LP(k,i) * LP(l,j);
               end
            end
            kof = sqrt((2*k-1)*(2*l-1))/2;      % including 0  term
            %     kof = sqrt((2*k+1)*(2*l+1))/2;
            
            CovChi(k,l) = CovChi(k,l)*SizeDomain* ...
               kof;
            CovChi(l,k) = CovChi(k,l);		   
         end
      end
      %       Compute the eigenvalues in order to transform Chi into Xi,
      %       ie independent variables
      options.disp = 0;
      %  Verifies CovMat*Phi =  Phi*Theta
      [Phi , Lambda] = eigs(CovChi,eye(size(CovChi)),OrderGaussExp,options);
      RF.Phi = Phi;              % eigenvectors of correlation/and
      % covariance matrix 
      RF.Eigs = diag(Lambda);  % eigenvalues of CorrMatrix order from
      
      RF.CovChi = CovChi ; 
      
   otherwise
      fprintf('\n   OSE still not available for 2D random fields\n');
      error(' ');
   end;
   
otherwise
   fprintf('\n\n NB : Valid discretization schemes are :\n');
   fprintf('          ''KL'',''EOLE'',''OSE''\n'); 
   error(' ');
end;

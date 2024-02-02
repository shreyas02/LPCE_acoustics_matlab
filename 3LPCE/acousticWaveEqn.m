%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                Acoustic equation solver                                 %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol, AcNormIndicator] = acousticWaveEqn(Sol, solver, acoustic, pmc, cnn, crd, elemType,...
                            nen, ndof, nElem, BCFluid, acousticBCNodes, ...
                            acousticBCValues)

% Quadrature rules for elements
if strcmp(elemType,'3Tri')
    gP = ...
   [1/3,  1/3
    4/3,  1/3
    1/3,  4/3] ;
    gW = ...
   [2/3,  2/3,  2/3] ;
 
    N(:,1) = 0.5.*(2.-gP(:,1)-gP(:,2)) ;
    N(:,2) = 0.5.*(gP(:,1)) ;
    N(:,3) = 0.5.*(gP(:,2)) ;
    
    Nx(:,1) = -0.5.*ones(3,1) ;
    Nx(:,2) =  0.5.*ones(3,1) ;
    Nx(:,3) =  zeros(3,1) ; 
    Ny(:,1) = -0.5.*ones(3,1) ;
    Ny(:,2) =  zeros(3,1) ;
    Ny(:,3) =  0.5.*ones(3,1) ;    
elseif strcmp(elemType,'4Quad')
    gP = ...
   [-5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01] ;
    gW = [1, 1, 1, 1 ] ;
    
    N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
    N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
    N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
    N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
    
    Nx(:,1) = -0.25.*(1-gP(:,2)) ;
    Nx(:,2) =  0.25.*(1-gP(:,2)) ;
    Nx(:,3) =  0.25.*(1+gP(:,2)) ;
    Nx(:,4) = -0.25.*(1+gP(:,2)) ;
    Ny(:,1) = -0.25.*(1-gP(:,1)) ;
    Ny(:,2) = -0.25.*(1+gP(:,1)) ;
    Ny(:,3) =  0.25.*(1+gP(:,1)) ;
    Ny(:,4) =  0.25.*(1-gP(:,1)) ;
end

Nx = Nx' ;
Ny = Ny' ;
nQuad = length(gW) ;

% Inverse mapping from global crd to structural crd as well as cnn
%crdMap = [[1:size(crdS,1)]' unique(BCStructure(:))];
%for i=1:size(cnnS,1)
%    for j=1:size(cnnS,2)
%        aa = find(crdMap(:,2) == cnnS(i,j));
%        cnnMap(i,j) = aa;
%    end
%end

pDashDot = Sol.pDashDot;
pDashDotPrev = Sol.pDashDotPrev;
pDashDDot = Sol.pDashDDot;
pDashDDotPrev = Sol.pDashDDotPrev;
pDash = Sol.pDash;
pDashPrev = Sol.pDashPrev;

% Form the local to global map
iif = zeros(nen^2*nElem,1); 
jjf = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElem) = double(cnn(:,i)); 
      jjf(index+1:index+nElem) = double(cnn(:,j));  
      index = index + nElem;
   end
end

% Satisfy boundary conditions
pDash(acousticBCNodes,1) = acousticBCValues ;

% Acoustic source term
acoustic.source = zeros(ndof,1);

% Interpolate for alpha values for Gen-alpha
pDashAlpha = pDashPrev + pmc.alpha.*(pDash - pDashPrev) ;
pDashDDotAlpha = pDashDDotPrev + pmc.alphaM.*(pDashDDot - pDashDDotPrev) ;
        
% Acoustic equation
xx = zeros(size(cnn));
yy = zeros(size(cnn));

% Localize the data to each element
for i=1:nen
   xx(:,i) = crd(cnn(:,i),1);
   yy(:,i) = crd(cnn(:,i),2);
   pD(:,i) = pDashAlpha(cnn(:,i),1);
end

% Form element matrix and assemble Galerkin terms
sA1 = zeros(nen^2*nElem,nQuad); 
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

for p = 1:nQuad  
    J = [xx*[Nx(:,p)], xx*[Ny(:,p)],...
         yy*[Nx(:,p)], yy*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
           
    volume = abs(volume);

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
               
    index = 0;
    for i = 1:nen
        for j = 1:nen    
            % Galerkin inertia term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = Mij.*(1/(acoustic.C*acoustic.C)) ;
            Mij = Mij.*volume;
            sA(index+1:index+nElem,p) = Mij;
            
            % Galerkin diffusion term
            Kij_1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Kij_2 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Kij_1 = Kij_1.*volume;
            Kij_2 = Kij_2.*volume;
            sK_1(index+1:index+nElem,p) = Kij_1;
            sK_2(index+1:index+nElem,p) = Kij_2;

            % Galerkin source term
            Aij_5 = gW(p)*(N(p,i)*N(p,j));
            Aij_5 = Aij_5.*volume ;
            sA5(index+1:index+nElem,p) = Aij_5;

            index = index + nElem;
        end
    end
end
% Summation of all quadrature data
sA = sum(sA,2);
sK_1 = sum(sK_1,2);
sK_2 = sum(sK_2,2);
sA5 = sum(sA5,2);

% Assemble the matrix
M = sparse(iif,jjf,sA,ndof,ndof);  
        
K1 = sparse(iif,jjf,sK_1,ndof,ndof);
K2 = sparse(iif,jjf,sK_2,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);

Src = [A5*acoustic.source(:)];
   
MAcoustic = (pmc.alphaM/(pmc.beta*pmc.alpha*(solver.dt)^2))*M;

% Left-hand side matrix
LHS = MAcoustic + K1 + K2;
        
clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
clear sA10 sA11 sA12 sA13 sA14 sA15 sA16
clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12 Aij_13 Aij_14 Aij_15 Aij_16
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16


% Right-hand side vector
%RHS = -(MAcoustic * pDashDDotAlpha(:)) - (K1+K2)*pDashAlpha(:);
RHS = -(MAcoustic * pDashAlpha(:)) + MAcoustic *(pDashPrev(:)) ...
    + (pmc.alphaM/(pmc.beta*solver.dt))* M * pDashDotPrev(:) ...
    - (1 - pmc.alphaM/(2*pmc.beta))* M * pDashDDotPrev(:) - (K1+K2)*pDashAlpha(:) ;
RHS = RHS + Src ;


% Solve the linear system
        
% Select the unknown nodal values
freeNodes = [acousticBCNodes] ;
freeNodes = setdiff(unique(BCFluid(:)),[freeNodes]);
        
result = pDashAlpha(:,1);
result = result(:);

Increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);
        
% Update the increments
result(freeNodes) = result(freeNodes) + Increment;

pDashAlpha(:,1) = result;

% Update the solution
pDash = pDashPrev + (1/pmc.alpha)*( pDashAlpha - pDashPrev ) ;
pDashDot = (pmc.gamma/(pmc.beta*solver.dt))* ( pDash - pDashPrev ...
     - solver.dt.*(1 - pmc.beta/pmc.gamma)*pDashDotPrev ...
     - solver.dt^2 * (0.5 - pmc.beta/pmc.gamma)*pDashDDotPrev ) ;
pDashDDot = (1/pmc.gamma).*( (pDashDot - pDashDotPrev)/solver.dt - (1-pmc.gamma)*pDashDDotPrev ) ;
        
AcNormIndicator =  norm(Increment)/(norm(result(freeNodes))+1e-12) ;
fprintf('Acoustic norm: %e, ', AcNormIndicator);
clear freeNodes1 freeNodes2 freeNodes3
clear result resultDot

% Map local structure to global data
Sol.pDashDot = pDashDot;
Sol.pDashDotPrev = pDashDotPrev;
Sol.pDashDDot = pDashDDot;
Sol.pDashDDotPrev = pDashDDotPrev;
Sol.pDash = pDash;
Sol.pDashPrev = pDashPrev;

end

% Petrov-Galerkin stabilization terms for Navier-Stokes equations
% Part 1: (tauM/rho)(rho*u.gradN).(rho*dU/dt) + 
%         (tauM/rho)(rho*u.gradN).(rho*u.gradU) +
%         (tauM/rho)(gradN).(rho*dU/dt)
%

function [LHS, RHS] = PetrovGalerkinTerms1(Sol, fluid, pmc, solver, iif, jjf, ...
                                          xxf, yyf, ux, uy, gW, N, Nx, Ny, ...
                                          nElem, nQuad, nen, ndof, elemType, ...
                                          LHS, RHS, alevelx, alevely)
                                      
                                      
sA1 = zeros(nen^2*nElem,nQuad); 
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

sA4 = zeros(nen^2*nElem,nQuad);
sA5 = zeros(nen^2*nElem,nQuad);
sA6 = zeros(nen^2*nElem,nQuad);      
sA7 = zeros(nen^2*nElem,nQuad); 
sA8 = zeros(nen^2*nElem,nQuad);
sA9 = zeros(nen^2*nElem,nQuad);
sA10 = zeros(nen^2*nElem,nQuad);
sA11 = zeros(nen^2*nElem,nQuad);
sA12 = zeros(nen^2*nElem,nQuad);
        
sA13 = zeros(nen^2*nElem,nQuad); 
sA14 = zeros(nen^2*nElem,nQuad);
sA15 = zeros(nen^2*nElem,nQuad);

for p = 1:nQuad  
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
    xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
    xsInv(:,4) = 1./xsInv(:,4) ;
    xsInv(:,1) = xsInv(:,4).* J(:,4) ;
    xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
    xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
    xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
    if (strcmp(elemType,'3Tri'))
        gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
        gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
        gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
    elseif (strcmp(elemType,'4Quad'))
        gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
        gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
        gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
    end
         
    volume = abs(volume);
    
    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
    locUX  = sum(repmat(N(p,:),nElem,1).*(ux - alevelx),2);
    locUY  = sum(repmat(N(p,:),nElem,1).*(uy - alevely),2); 
                
    VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
    VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

    VGV	= VG1 .* locUX + VG2 .* locUY ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
            
    tauM = (2/solver.dt)^2 + ((fluid.visc./fluid.dens).^2).* GG * 36 + abs(VGV);    
    tauM = tauM.^(0.5);
    tauM = 1./tauM;
            
    index = 0;
    for i = 1:nen
        for j = 1:nen
            % 1st term
            Aij_1 = gW(p)*(locUX.*DNDx(:,i)*N(p,j));
            Aij_2 = gW(p)*(locUY.*DNDy(:,i)*N(p,j));
            Aij_1 = Aij_1.*volume.*tauM.*fluid.dens;
            Aij_2 = Aij_2.*volume.*tauM.*fluid.dens;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
                    
            % 2nd term
            Aij_4 = gW(p)*locUX.*locUX.*(DNDx(:,i).*DNDx(:,j));
            Aij_5 = gW(p)*locUX.*locUY.*(DNDx(:,i).*DNDy(:,j));
            Aij_7 = gW(p)*locUY.*locUX.*(DNDy(:,i).*DNDx(:,j));
            Aij_8 = gW(p)*locUY.*locUY.*(DNDy(:,i).*DNDy(:,j)); 
            Aij_4 = Aij_4.*volume.*tauM.*fluid.dens;
            Aij_5 = Aij_5.*volume.*tauM.*fluid.dens;
            Aij_7 = Aij_7.*volume.*tauM.*fluid.dens;
            Aij_8 = Aij_8.*volume.*tauM.*fluid.dens;
            sA4(index+1:index+nElem,p) = Aij_4;
            sA5(index+1:index+nElem,p) = Aij_5;
            sA7(index+1:index+nElem,p) = Aij_7;
            sA8(index+1:index+nElem,p) = Aij_8;

            % 3rd term
            Aij_13 = gW(p)*(DNDx(:,i).*N(p,j));
            Aij_14 = gW(p)*(DNDy(:,i).*N(p,j));
            Aij_13 = Aij_13.*volume.*tauM;
            Aij_14 = Aij_14.*volume.*tauM; 
            sA13(index+1:index+nElem,p) = Aij_13;
            sA14(index+1:index+nElem,p) = Aij_14;
      
            index = index + nElem;
        end
    end
end
% Summation of all quadrature data
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);       
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA7 = sum(sA7,2);
sA8 = sum(sA8,2);
sA13 = sum(sA13,2);
sA14 = sum(sA14,2);
        
% Assemble the matrix        
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A7 = sparse(iif,jjf,sA7,ndof,ndof);
A8 = sparse(iif,jjf,sA8,ndof,ndof);
A13 = sparse(iif,jjf,sA13,ndof,ndof);
A14 = sparse(iif,jjf,sA14,ndof,ndof);
        
ZeroF = sparse(ndof,ndof);

stabK1 = [A1+A2 ZeroF ZeroF;...
          ZeroF A1+A2 ZeroF;...
          ZeroF ZeroF ZeroF];
stabK2 = [A4+A5+A7+A8 ZeroF ZeroF;...
          ZeroF A4+A5+A7+A8 ZeroF;...
          ZeroF ZeroF ZeroF];
stabG1 = [ZeroF ZeroF ZeroF;...
          ZeroF ZeroF ZeroF;...
          A13 A14 ZeroF];

stabK1Flow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabK1;
stabG1Flow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabG1;

% Left-hand side matrix
LHS = LHS + (stabK1Flow + stabK2 + stabG1Flow);
        
clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 sA13 sA14 sA15
clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
clear Aij_13 Aij_14 Aij_15
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 xsInv

uAlpha = Sol.uAlpha(:,:,1) ;
uAlpha = [uAlpha(:); zeros(ndof,1)];
uDotAlpha = Sol.uDotAlpha(:,:,1) ;
uDotAlpha = [uDotAlpha(:); zeros(ndof,1)];
if (size(Sol.p,1) == 1 & size(Sol.p,2)>1)
   Sol.p = [Sol.p]';
end

p1 = [zeros(ndof*2,1);Sol.p];
gravVec = [ones(2*ndof,1); zeros(ndof,1)];
        
% Right-hand side vector
RHS = RHS -(stabK1 + stabG1) * uDotAlpha(:) - (stabK2) * uAlpha(:) ;                                   
                                      
end
% Petrov-Galerkin stabilization terms for Navier-Stokes equations
% Part 3: (tauM/rho)(gradN).(rho*u.(gradU)) + 
%         (tauM/rho)(gradN).(gradP) +
%         (tauM/rho)(gradN).(-rho.gravVec)
%

function [LHS, RHS] = PetrovGalerkinTerms3(Sol, fluid, pmc, solver, iif, jjf, ...
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
            Aij_1 = gW(p)*locUX.*(DNDx(:,i).*DNDx(:,j));
            Aij_2 = gW(p)*locUY.*(DNDx(:,i).*DNDy(:,j));
            Aij_4 = gW(p)*locUX.*(DNDy(:,i).*DNDx(:,j));
            Aij_5 = gW(p)*locUY.*(DNDy(:,i).*DNDy(:,j));
            Aij_1 = Aij_1.*volume.*tauM;
            Aij_2 = Aij_2.*volume.*tauM;
            Aij_4 = Aij_4.*volume.*tauM;
            Aij_5 = Aij_5.*volume.*tauM;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
            sA4(index+1:index+nElem,p) = Aij_4;
            sA5(index+1:index+nElem,p) = Aij_5;

            % 2nd term
            Aij_10 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Aij_11 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Aij_10 = Aij_10.*volume.*tauM.*(1./fluid.dens);
            Aij_11 = Aij_11.*volume.*tauM.*(1./fluid.dens);
            sA10(index+1:index+nElem,p) = Aij_10;
            sA11(index+1:index+nElem,p) = Aij_11;
                    
            % 3rd term
            Aij_13 = gW(p)*(DNDx(:,i)*N(p,j));
            Aij_14 = gW(p)*(DNDy(:,i)*N(p,j));
            Aij_13 = Aij_13.*volume.*tauM ;
            Aij_14 = Aij_14.*volume.*tauM ;
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
sA10 = sum(sA10,2);
sA11 = sum(sA11,2);
sA13 = sum(sA13,2);
sA14 = sum(sA14,2);
        
% Assemble the matrix        
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A10 = sparse(iif,jjf,sA10,ndof,ndof);
A11 = sparse(iif,jjf,sA11,ndof,ndof);
A13 = sparse(iif,jjf,sA13,ndof,ndof);
A14 = sparse(iif,jjf,sA14,ndof,ndof);
       
ZeroF = sparse(ndof,ndof);
 
stabG2 = [ZeroF ZeroF ZeroF;...
          ZeroF ZeroF ZeroF;...
          A1+A2 A4+A5 ZeroF];
stabC1 = [ZeroF ZeroF ZeroF;...
          ZeroF ZeroF ZeroF;...
          ZeroF ZeroF A10+A11];
Src4 = [ZeroF ZeroF ZeroF;...
        ZeroF ZeroF ZeroF;...
        A13.*fluid.gravFrc(1) A14.*fluid.gravFrc(2) ZeroF] ;

% Left-hand side matrix
LHS = LHS + (stabG2 + stabC1);
        
clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 
clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12
        
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
RHS = RHS - stabG2 * uAlpha(:) ;
RHS = RHS - stabC1 * p1 + Src4 * gravVec;                                    
                                  
end
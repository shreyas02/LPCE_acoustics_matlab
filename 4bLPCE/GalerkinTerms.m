% Galerkin terms for Navier-Stokes equations

function [LHS, RHS] = GalerkinTerms(Sol, fluid, pmc, solver, iif, jjf, xxf, yyf, ux, uy, ...
                                    gW, N, Nx, Ny, nElem, nQuad, nen, ndof, alevelx,...
                                    alevely)

sA = [] ;
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
sA16 = zeros(nen^2*nElem,nQuad);

for p = 1:nQuad  
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
           
    volume = abs(volume);
    
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
       exit
    end

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
    locUX  = sum(repmat(N(p,:),nElem,1).*(ux-alevelx),2);
    locUY  = sum(repmat(N(p,:),nElem,1).*(uy-alevely),2);   
            
    index = 0;
    for i = 1:nen
        for j = 1:nen
            % Galerkin inertia term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = Mij.*fluid.dens ;
            Mij = Mij.*volume;
            sA(index+1:index+nElem,p) = Mij;
                    
            % Galerkin convection term
            Aij_1 = gW(p)*N(p,i)*(locUX.*DNDx(:,j));
            Aij_2 = gW(p)*N(p,i)*(locUY.*DNDy(:,j));
            Aij_1 = Aij_1.*volume.*fluid.dens;
            Aij_2 = Aij_2.*volume.*fluid.dens;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
                    
            % Galerkin viscous diffusion term
            Aij_4 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Aij_5 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Aij_7 = gW(p)*(DNDx(:,i).*DNDy(:,j));
            Aij_4 = Aij_4.*volume.*fluid.visc;
            Aij_5 = Aij_5.*volume.*fluid.visc;
            Aij_7 = Aij_7.*volume.*fluid.visc;
            sA4(index+1:index+nElem,p) = Aij_4;
            sA5(index+1:index+nElem,p) = Aij_5;
            sA7(index+1:index+nElem,p) = Aij_7;
                    
            % Galerkin pressure term (gradP)
            Aij_10 = gW(p)*(DNDx(:,i).*N(p,j));
            Aij_11 = gW(p)*(DNDy(:,i).*N(p,j));
            Aij_10 = Aij_10.*volume;
            Aij_11 = Aij_11.*volume; 
            sA10(index+1:index+nElem,p) = Aij_10;
            sA11(index+1:index+nElem,p) = Aij_11;
                    
            % Galerkin source term (gravity force)
            Aij_13 = gW(p)*(N(p,i)*N(p,j));
            Aij_13 = Aij_13.*fluid.dens.*volume ;
            sA13(index+1:index+nElem,p) = Aij_13;
                    
            % Galerkin continuity term
            Aij_14 = gW(p)*(N(p,i)).*(DNDx(:,j)) ;
            Aij_15 = gW(p)*(N(p,i)).*(DNDy(:,j)) ;
            Aij_14 = Aij_14.*volume ;
            Aij_15 = Aij_15.*volume ;
            sA14(index+1:index+nElem,p) = Aij_14;
            sA15(index+1:index+nElem,p) = Aij_15;
            
            index = index + nElem;
        end
    end
end
% Summation of all quadrature data
sA = sum(sA,2);
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA7 = sum(sA7,2);
sA10 = sum(sA10,2);
sA11 = sum(sA11,2);
sA13 = sum(sA13,2);
sA14 = sum(sA14,2);
sA15 = sum(sA15,2);
 
% Assemble the matrix
Mf = sparse(iif,jjf,sA,ndof,ndof);  
ZeroF = sparse(ndof,ndof);
        
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A7 = sparse(iif,jjf,sA7,ndof,ndof);
A10 = sparse(iif,jjf,sA10,ndof,ndof);
A11 = sparse(iif,jjf,sA11,ndof,ndof);
A13 = sparse(iif,jjf,sA13,ndof,ndof);
A14 = sparse(iif,jjf,sA14,ndof,ndof);
A15 = sparse(iif,jjf,sA15,ndof,ndof);

Mf = [Mf ZeroF ZeroF;...
      ZeroF Mf ZeroF;...
      ZeroF ZeroF ZeroF];
Conv = [A1+A2 ZeroF ZeroF;...
        ZeroF A1+A2 ZeroF;...
        ZeroF ZeroF ZeroF];
Kf = [2*A4+A5 A7' ZeroF;...
      A7 A4+2*A5  ZeroF;...
      ZeroF ZeroF ZeroF];
Fp = [ZeroF ZeroF A10;...
      ZeroF ZeroF A11;...
      ZeroF ZeroF ZeroF];
GMassCons = [ZeroF ZeroF ZeroF;...
             ZeroF ZeroF ZeroF;...
             A14 A15 ZeroF];
Src = [A13.*fluid.gravFrc(1) ZeroF ZeroF;...  
       ZeroF A13.*fluid.gravFrc(2) ZeroF;...;...
       ZeroF ZeroF ZeroF];
   
Mflow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Mf;

% Left-hand side matrix
LHS = Mflow + Conv  + Kf - Fp + GMassCons;
        
clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
clear sA10 sA11 sA12 sA13 sA14 sA15 sA16
clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12 Aij_13 Aij_14 Aij_15 Aij_16
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16

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
RHS = -(Mf * uDotAlpha(:)) - (Conv + Kf + GMassCons)* uAlpha(:) ;
RHS = RHS + (Fp * p1) + (Src)*gravVec ;

end
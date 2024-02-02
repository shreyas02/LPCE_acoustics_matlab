function [Sol, lpceNormIndicator] = lpce(Sol, solver, acoustic, pmc, cnn, crd, elemType,...
                            nen, ndof, BCFluid , fluid,nElem)

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


% Form the local to global map
ii = zeros(nen^2*nElem,1); 
jj = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      ii(index+1:index+nElem) = double(cnn(:,i)); 
      jj(index+1:index+nElem) = double(cnn(:,j));  
      index = index + nElem;
   end
end

% Satisfy boundary conditions
pa(acoustic.Dirichletp) =  acoustic.Dirichletpa;
ua(acoustic.Dirichletu) = acoustic.Dirichletua;
va(acoustic.Dirichletv) = acoustic.Dirichletva;
rho(acoustic.Dirichletrho) = acoustic.Dirichletrhoa;

%incorporating the incompressible flow data at n + alpha
ic.u = Sol.u(:,1);
ic.v = Sol.u(: ,2);
ic.p = Sol.p;

% Interpolate for alpha values for Gen-alpha

% Initialization for Element-level Calculations 
xx = zeros(size(cnn));
yy = zeros(size(cnn));

for i=1:nen
   xx(:,i) = crd(cnn(:,i),1);
   yy(:,i) = crd(cnn(:,i),2);
end

%Area of the elements
area = (xx(:,1) - xx(:,3)).*(yy(:,2) - yy(:,4)) - (xx(:,2) - xx(:,4)).*(yy(:,1) - yy(:,3));
area = 0.5.*abs(area);
%element length 
h = sqrt(area);

%Element matrix  incomopressible flow data
for i=1:nen
   ubar(:,i) = ic.u(cnn(:,i),1) ;
   vbar(:,i) = ic.v(cnn(:,i),1) ;
   pbar(:,i) = ic.p(cnn(:,i),1) ;
end

%Element-level Evaluation of Matrices

% Gauss quadrature loop
for p = 1:nQuad  
    
    % Jacobian evaluation and its absolute value
    J = [xx*[Nx(:,p)], yy*[Nx(:,p)],...
         xx*[Ny(:,p)], yy*[Ny(:,p)]];
    absJac =abs( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
    
    % Evaluation of gradients of shape functions 
    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,2))*Ny(:,p)')./repmat(absJac,1,nen);
    DNDy = ((-J(:,3))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(absJac,1,nen);
    
    uic  = sum(repmat(N(p,:),nElem,1).*(ubar),2);
    vic  = sum(repmat(N(p,:),nElem,1).*(vbar),2);
    pic = sum(repmat(N(p,:),nElem,1).*(pbar),2);   

    dudxic  = sum(DNDx.*(ubar),2);
    dudyic  = sum(DNDy.*(ubar),2);
    dvdxic  = sum(DNDx.*(vbar),2);
    dvdyic  = sum(DNDy.*(vbar),2);
    dpdxic = sum(DNDx.*(pbar),2);
    dpdyic = sum(DNDy.*(pbar),2);

    %stabalistation factor 
    beta = 1;
    tau = (beta.*h)./(2.*sqrt(uic.^2+vic.^2));

    % Evaluation of w tilda     
    wTilda = DNDx.*uic + DNDy.*vic;
    wTilda = tau.*wTilda;
            
    index = 0;
    for i = 1:nen
        for j = 1:nen
            
            %N_xN Galerkin continuity term
            Xij_1 = gW(p)*(N(p,i)).*(DNDx(:,j))  + gW(p)*(wTilda(:,i)).*(DNDx(:,j));
            Xij_1 = Xij_1.*absJac;
            xK1(index+1:index+nElem,p) = Xij_1;

            Xij_2 = gW(p)*(N(p,i)).*(DNDx(:,j).*uic)  + gW(p)*(wTilda(:,i)).*(DNDx(:,j).*uic);
            Xij_2 = Xij_2.*absJac;
            xK2(index+1:index+nElem,p) = Xij_2;

            Xij_3 = gW(p)*(N(p,i)).*(DNDx(:,j).*vic)  + gW(p)*(wTilda(:,i)).*(DNDx(:,j).*vic);
            Xij_3 = Xij_3.*absJac;
            xK3(index+1:index+nElem,p) = Xij_3;

            Xij_4 = gW(p)*(N(p,i)).*(DNDx(:,j).*pic) + gW(p)*(wTilda(:,i)).*(DNDx(:,j).*pic);
            Xij_4 = Xij_4.*absJac;
            xK4(index+1:index+nElem,p) = Xij_4;

            %N_yN Galerkin continuity term
            Yij_1 = gW(p)*(N(p,i)).*(DNDy(:,j)) + gW(p)*(wTilda(:,i)).*(DNDy(:,j));
            Yij_1 = Yij_1.*absJac;
            yK1(index+1:index+nElem,p) = Yij_1;

            Yij_2 = gW(p)*(N(p,i)).*(DNDy(:,j).*uic) + gW(p)*(wTilda(:,i)).*(DNDy(:,j).*uic);
            Yij_2 = Yij_2.*absJac;
            yK2(index+1:index+nElem,p) = Yij_2;

            Yij_3 = gW(p)*(N(p,i)).*(DNDy(:,j).*vic) + gW(p)*(wTilda(:,i)).*(DNDy(:,j).*vic);
            Yij_3 = Yij_3.*absJac;
            yK3(index+1:index+nElem,p) = Yij_3;

            Yij_4 = gW(p)*(N(p,i)).*(DNDy(:,j).*pic) + gW(p)*(wTilda(:,i)).*(DNDy(:,j).*pic);
            Yij_4 = Yij_4.*absJac;
            yK4(index+1:index+nElem,p) = Yij_4;

            %Galerkin mass matrix term 
            Mij_1 = gW(p)*(N(p,i)*N(p,j)) + gW(p)*(wTilda(:,i)).*N(p,j);
            Mij_1 = Mij_1.*absJac;
            mK1(index+1:index+nElem,p) = Mij_1;

            Mij_2 = gW(p)*(N(p,i)*N(p,j).*dudxic) + gW(p)*(wTilda(:,i)).*N(p,j).*dudxic;
            Mij_2 = Mij_2.*absJac;
            mK2(index+1:index+nElem,p) = Mij_2;

            Mij_3 = gW(p)*(N(p,i)*N(p,j).*dvdyic) + gW(p)*(wTilda(:,i)).*N(p,j).*dvdyic;
            Mij_3 = Mij_3.*absJac;
            mK3(index+1:index+nElem,p) = Mij_3;

            Mij_4 = gW(p)*(N(p,i)*N(p,j).*dvdxic) + gW(p)*(wTilda(:,i)).*N(p,j).*dvdxic;
            Mij_4 = Mij_4.*absJac;
            mK4(index+1:index+nElem,p) = Mij_4;

            Mij_5 = gW(p)*(N(p,i)*N(p,j).*dudyic) + gW(p)*(wTilda(:,i)).*N(p,j).*dudyic;
            Mij_5 = Mij_5.*absJac;
            mK5(index+1:index+nElem,p) = Mij_5;

            Mij_6 = gW(p)*(N(p,i)*N(p,j).*dpdxic) + gW(p)*(wTilda(:,i)).*N(p,j).*dpdxic;
            Mij_6 = Mij_6.*absJac;
            mK6(index+1:index+nElem,p) = Mij_6;

            Mij_7 = gW(p)*(N(p,i)*N(p,j).*dpdyic) + gW(p)*(wTilda(:,i)).*N(p,j).*dpdyic;
            Mij_7 = Mij_7.*absJac;
            mK7(index+1:index+nElem,p) = Mij_7;

            index = index + nElem;
        end
    end
end

% Summation of All Quadrature Data for Numerical Integration
mK1=sum(mK1,2);
mK2=sum(mK2,2);
mK3=sum(mK3,2);
mK4=sum(mK4,2);
mK5=sum(mK5,2);
mK6=sum(mK6,2);
mK7=sum(mK7,2);

xK1=sum(xK1,2);
xK2=sum(xK2,2);
xK3=sum(xK3,2);
xK4=sum(xK4,2);

yK1=sum(yK1,2);
yK2=sum(yK2,2);
yK3=sum(yK3,2);
yK4=sum(yK4,2);

% Assembly of Local Matrices to Global Matrix

M1 = sparse(ii,jj,mK1,ndof,ndof); 
M2 = sparse(ii,jj,mK2,ndof,ndof); 
M3 = sparse(ii,jj,mK3,ndof,ndof); 
M4 = sparse(ii,jj,mK4,ndof,ndof); 
M5 = sparse(ii,jj,mK5,ndof,ndof); 
M6 = sparse(ii,jj,mK6,ndof,ndof); 
M7 = sparse(ii,jj,mK7,ndof,ndof); 

X1 = sparse(ii,jj,xK1,ndof,ndof); 
X2 = sparse(ii,jj,xK2,ndof,ndof);
X3 = sparse(ii,jj,xK3,ndof,ndof); 
X4 = sparse(ii,jj,xK4,ndof,ndof); 

Y1 = sparse(ii,jj,yK1,ndof,ndof); 
Y2 = sparse(ii,jj,yK2,ndof,ndof); 
Y3 = sparse(ii,jj,yK3,ndof,ndof); 
Y4 = sparse(ii,jj,yK4,ndof,ndof); 

Zero2  = sparse(ndof,ndof);
Zero1 = sparse(ndof,1);

% Select the unknown nodal values
freeNodesu = [acoustic.Dirichletu] ;
freeNodesu = setdiff(unique(BCFluid(:)),[freeNodesu]);

freeNodesv = [acoustic.Dirichletv] ;
freeNodesv = setdiff(unique(BCFluid(:)),[freeNodesv]);

freeNodesp = [acoustic.Dirichletp] ;
freeNodesp = setdiff(unique(BCFluid(:)),[freeNodesp]);

freeNodesrho = [acoustic.Dirichletrho] ;
freeNodesrho = setdiff(unique(BCFluid(:)),[freeNodesrho]);

freeNodes = [freeNodesrho ;freeNodesu + ndof; freeNodesv + 2*ndof; freeNodesp + 3*ndof];

% Initialising the forcing function 

Q = [Zero1 ; Zero1 ; Zero1 ; M1*Sol.q];

% Construction of the final matrices  
A = [M1, Zero2 , Zero2 , Zero2 ;...
        Zero2, M1 , Zero2, Zero2 ;...
        Zero2, Zero2, M1, Zero2;...
        Zero2, Zero2, Zero2, M1];

B = [X2 + Y3 , fluid.dens.*X1 , fluid.dens.*Y1, Zero2;...
        Zero2, X2 + M2, X3 + M4, (1/fluid.dens).*X1;...
        Zero2, Y2 + M5, Y3 + M3, (1/fluid.dens).*Y1;...
        Zero2, fluid.gamma.*X4 + M6, fluid.gamma.*Y4 + M7, X2 + Y3];

%previous time step
oldResult = [Sol.rhoPrev ; Sol.uaPrev; Sol.vaPrev; Sol.paPrev];

%Defining coefficients 
coeffA = pmc.alphaM/(pmc.alpha*pmc.gamma*solver.dt);
coeffB = (1 - pmc.alphaM/pmc.gamma);

LHS = coeffA.*A + B;
RHS = Q + coeffA.*A*oldResult - coeffB*A*acoustic.v;

Result(freeNodes) = LHS(freeNodes,freeNodes)\RHS(freeNodes);

%transforming result - 
alphaResult = Result';
Result = reshape(Result , ndof , 4);
rho = Result(: , 1);
ua = Result(: , 2);
va = Result(: , 3);
pa = Result(: , 4);

% N+1 result from n + alpha
rho = (1/(pmc.alpha)).*(rho - (1-pmc.alpha).*Sol.rhoPrev);
ua = (1/(pmc.alpha)).*(ua - (1-pmc.alpha).*Sol.uaPrev);
va = (1/(pmc.alpha)).*(va - (1-pmc.alpha).*Sol.vaPrev);
pa = (1/(pmc.alpha)).*(pa - (1-pmc.alpha).*Sol.paPrev);


% Update acoustic.v 
acoustic.vAlpham = coeffB.*acoustic.v + coeffA.*(alphaResult - oldResult);
acoustic.v = (acoustic.vAlpham - (1-pmc.alphaM).*acoustic.v)./pmc.alphaM;

% Norm error generated
lpceNormIndicator = norm(Sol.paPrev - pa)/(norm(pa) + 1e-12);

%Introducing NRBC boundary conditions 
for k = 1 : ndof
    if(Sol.r(k) > 100)
        xb = (Sol.r(k) - 100)./50;
        c1 = 0.05;
        c2 = 20;
        damp = (1 - c1*xb^2) * (1 - ((1 - exp(c2 * xb^2))/(1 - exp(c2))));
        rho(k) = damp * rho(k);
        ua(k) = damp * ua(k);
        va(k) = damp * va(k);
        pa(k) = damp * pa(k);
   end
end

% Map local structure to global data 
Sol.rho = rho; 
Sol.ua = ua; 
Sol.va = va; 
Sol.pa = pa; 

%print error 
fprintf('LPCE norm: %e, ', lpceNormIndicator);

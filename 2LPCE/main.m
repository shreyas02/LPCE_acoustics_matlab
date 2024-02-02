%%%%% 1LPCE SOLVER %%%%% 
clc
clear all
wrkDir = './' ;
problemString = 'acoustic' ;
elemType = '4Quad' ;
problemType = '2D' ;

% Coordinate, connectivity and boundary data
dataStr = strcat(wrkDir,'Mesh.mat') ;
load(dataStr);
ndof = size(crd,1) ;

% Nonlinear iteration data:
solver.nLIterMin = 2 ;
solver.nLIterMax = 1;
solver.dt = 1/3740 ;
solver.maxSteps = 2750; %2750
solver.rhoinftyF = 1.0 ;
solver.nLTol = 5e-4 ;
solver.outFreq = [1, 1] ; % After time step solver.outFreq(1), outFreq will be solver.outFreq(2)
solver.intgOutFreq = 1 ;
restartFlag = 0 ;

% Fluid properties:
Re = 100 ;
fluid.dens = 1.2 ;
fluid.visc = fluid.dens*2*1/Re ;
fluid.gravFrc = [0, 0];
r_0 = [0, 0];
fluid.gamma = 1.4;

% Acoustic properties:
acoustic.C = 331;
acoustic.M = 0.1;

% Initial boundary conditions:
fluid.vel0 = [0, 0];
fluid.pres0 = 0.0 ;

% Boundary conditions:
fluid.Bc_CylN = size(unique(BCCyl),1) ;
fluid.Bc_CylNs = unique(BCCyl);
fluid.Bc_CylV = [0 0];

fluid.Bc_OutercylN = size(unique(BCOutercyl),1) ;
fluid.Bc_OutercylNs = unique(BCOutercyl);
fluid.Bc_OutercylV = [0 0];

fluid.radialN = size(unique(BCradial),1) ;
fluid.radialNs = unique(BCradial);
fluid.Bc_CylV = [0 0];

fluid.Bc_WallN = [] ;
fluid.Bc_WallNs = [];
fluid.Bc_WallV = [0 0];

fluid.DirichletU = [];
fluid.DirichletUval = [];
fluid.DirichletV = [];
fluid.DirichletVval = [];

%%%%%Acoustic Boundary conditions 
acoustic.Dirichletu= [fluid.Bc_CylNs ; fluid.Bc_OutercylNs];
acoustic.Dirichletv= [fluid.Bc_CylNs ; fluid.Bc_OutercylNs];
acoustic.Dirichletrho=  [fluid.Bc_OutercylNs];
acoustic.Dirichletp=   [fluid.Bc_OutercylNs];
acoustic.Dirichletua = [zeros(size(acoustic.Dirichletu,1),1)];
acoustic.Dirichletva = [zeros(size(acoustic.Dirichletv,1),1)];
acoustic.Dirichletpa = [zeros(size(acoustic.Dirichletp,1),1)];
acoustic.Dirichletrhoa = [zeros(size(acoustic.Dirichletrho,1),1)];

pl = zeros(fluid.radialN,3);

%BCInner = [unique(BCCylinder(:)); unique(BCStructure(:))]; 
%meshBCNodes = [BCBound; BCInner];

oisdNodes = [] ;
fluidBCNodesU = [fluid.DirichletU; oisdNodes]; % Wall boundary condition is satisfied for oisdNodes
fluidBCNodesV = [fluid.DirichletV; oisdNodes];
fluidBCNodesP = [];
Force = zeros(size(oisdNodes,1),2);
Traction = [0 0];
Sol.Traction = [0 0];
Sol.TractionPrev = [0 0] ;

% Gen-alpha time integration parameters
pmc.alphaM = 0.5*(3-solver.rhoinftyF)/(1+solver.rhoinftyF) ;
pmc.alpha = 1/(1+solver.rhoinftyF) ;
pmc.gamma = 0.5 + pmc.alphaM - pmc.alpha ;
pmc.beta = 0.25*(1 + pmc.alphaM - pmc.alpha)^2 ;

% Initialize other variables
nodeId = crd(:,1) ;
crd = crd(:,2:4) ;
cnnF = BCFluid ;
nenF = size(cnnF,2) ;
nElemF = size(cnnF,1) ;
cnn = [cnnF]; 
ndof = size(crd,1);
nen = nenF ;

% Fluid variables
Sol.u = zeros(ndof,2,1);
Sol.uDot = zeros(ndof,2,1);
Sol.u(:,1,:) = fluid.vel0(1) ;
Sol.u(:,2,:) = fluid.vel0(2) ;
Sol.uAlpha = zeros(ndof,2,1) ;
Sol.uDotAlpha = zeros(ndof,2,1) ;
Sol.uPrev = Sol.u ;
Sol.uDotPrev = Sol.uDot ;
Sol.p = fluid.pres0*ones(ndof,1);
Sol.pPrev = Sol.p ;
Sol.rho = zeros(ndof,1);

% Acoustic variables
Sol.rho = zeros(ndof,1);
Sol.va = zeros(ndof,1);
Sol.ua = zeros(ndof,1);
Sol.pa = zeros(ndof,1);

Sol.rhoPrev = zeros(ndof,1);
Sol.vaPrev = zeros(ndof,1);
Sol.uaPrev = zeros(ndof,1);
Sol.paPrev = zeros(ndof,1);

%generalised alpha acoustic
acoustic.v = zeros(4*ndof, 1);

% exact solution prev
exact.vrPrev = zeros(ndof,1);
exact.uPrev = zeros(ndof,1);
exact.vPrev = zeros(ndof,1);
exact.pPrev = zeros(ndof,1);
exact.DPDTPrev = zeros(ndof,1);

% ALE mesh variables
Sol.aleDisp = zeros(ndof,2);
Sol.aleDispPrev = zeros(ndof,2);
Sol.aleDispPrevPrev = zeros(ndof,2);
Sol.aleVel = zeros(ndof,2);
Sol.aleVelPrev = zeros(ndof,2);

crdNew = crd;

% Load previous time step (for restart)
if restartFlag == 1
    load('TS_restart.mat');
    rsTsId = 2200 ;
    solver.maxSteps = 3000 ;
else
    rsTsId = 1;
end

%Exact Solution at t = 0 ; 
    
time = 0.0 ; 
exact.r = sqrt(crd(:,1).^2  + crd(:,2).^2);
Sol.r = exact.r; 
exact.u = acoustic.M * acoustic.C;
exact.omega = 0.2.* acoustic.C .* (1/0.5) .* 2*pi ; %Strouhal number = 0.2

exact.vrPrev = exact.u.*((0.5^2)./(exact.r.^2)).*sin(exact.omega * time);
exact.uPrev = exact.vrPrev.*(crd(:,1)./exact.r);
exact.vPrev = exact.vrPrev.*(crd(:,2)./exact.r);
exact.dpdtn = -fluid.dens.*exact.omega^2.*0.5^2.*exact.u.*(1./(exact.r)).*sin(exact.omega*time) - ...
                                    fluid.dens.*exact.omega.*0.5^4.*exact.u^2.*(1./(exact.r.^4)).*sin(exact.omega*time).*cos(exact.omega*time);
exact.dpdrn = -fluid.dens.*exact.omega.*0.5^2.*exact.u.*(1./(exact.r.^2)).*cos(exact.omega*time) +...
                                    2.*fluid.dens.*0.5^4.*exact.u^2.*(1./(exact.r.^5)).*sin(exact.omega*time).*sin(exact.omega*time);
exact.dpdxn = exact.dpdrn .*(crd(:,1)./exact.r);
exact.dpdyn = exact.dpdrn .*(crd(:,2)./exact.r);
exact.DPDTPrev = exact.dpdtn + exact.dpdxn.*exact.uPrev + exact.dpdyn.*exact.vPrev;
exact.pPrev =  ((fluid.dens.*acoustic.C^2)/fluid.gamma).*ones(ndof,1) + fluid.dens.*exact.omega.*0.5^2.*exact.u.*(1./exact.r).*cos(exact.omega*time) -...
                            fluid.dens.*0.5^4.*exact.u^2.*(1./(2*exact.r.^4)).*sin(exact.omega*time).*sin(exact.omega*time);

% Time loop starts here
filename1 = sprintf('%s/%s_%d.oisd',wrkDir,problemString, rsTsId);
fileId1 = fopen(filename1,'w');
filename2 = sprintf('%s/%s_%d.othd',wrkDir,problemString, rsTsId);
fileId2 = fopen(filename2,'w');
time = 0.0 ;

for timeStep = rsTsId:solver.maxSteps
    if (mod(timeStep,100) == 0 )
        save('TS_restart.mat');
    end

    fprintf('Time step:%d\n',timeStep);
    time = timeStep.*solver.dt ;
    timePrev = (timeStep-1).*solver.dt ;
    
    % Predict the solution
    Sol.u = Sol.uPrev ;
    Sol.uDot = (pmc.gamma - 1)/pmc.gamma * Sol.uDotPrev ;
    Sol.p = Sol.pPrev ;
    Sol.aleDisp = Sol.aleDispPrev ;
    Sol.aleDispPrev = Sol.aleDispPrevPrev ;
    Sol.aleVel = Sol.aleVelPrev ;

    % Predict the solution Acoustic 
    Sol.ua = Sol.uaPrev;
    Sol.va = Sol.vaPrev;
    Sol.pa = Sol.paPrev;
    Sol.rho = Sol.rhoPrev;

    %Exact Solution 
    exact.r = sqrt(crd(:,1).^2  + crd(:,2).^2);
    Sol.r = exact.r; 
    exact.u = acoustic.M * acoustic.C;
    exact.omega = 0.2.* acoustic.C .* (1/0.5) .* 2*pi ; %Strouhal number = 0.2

    %at time n + 1 
    exact.vrn = exact.u.*((0.5^2)./(exact.r.^2)).*sin(exact.omega * time);
    exact.un = exact.vrn.*(crd(:,1)./exact.r);
    exact.vn = exact.vrn.*(crd(:,2)./exact.r);
    exact.dpdtn = -fluid.dens.*exact.omega^2.*0.5^2.*exact.u.*(1./(exact.r)).*sin(exact.omega*time) - ...
                                    fluid.dens.*exact.omega.*0.5^4.*exact.u^2.*(1./(exact.r.^4)).*sin(exact.omega*time).*cos(exact.omega*time);
    exact.dpdrn = -fluid.dens.*exact.omega.*0.5^2.*exact.u.*(1./(exact.r.^2)).*cos(exact.omega*time) +...
                                    2.*fluid.dens.*0.5^4.*exact.u^2.*(1./(exact.r.^5)).*sin(exact.omega*time).*sin(exact.omega*time);
    exact.dpdxn = exact.dpdrn .*(crd(:,1)./exact.r);
    exact.dpdyn = exact.dpdrn .*(crd(:,2)./exact.r);
    exact.DPDTn = exact.dpdtn + exact.dpdxn.*exact.un + exact.dpdyn.*exact.vn;
    exact.pn =  ((fluid.dens.*acoustic.C^2)/fluid.gamma).*ones(ndof,1) + fluid.dens.*exact.omega.*0.5^2.*exact.u.*(1./exact.r).*cos(exact.omega*time) -...
                            fluid.dens.*0.5^4.*exact.u^2.*(1./(2*exact.r.^4)).*sin(exact.omega*time).*sin(exact.omega*time);

    %at time n + alpha 
    exact.vra = pmc.alpha.*exact.vrn + (1-pmc.alpha).*exact.vrPrev;
    exact.ua = pmc.alpha.*exact.un + (1-pmc.alpha).*exact.uPrev;
    exact.va = pmc.alpha.*exact.vn + (1-pmc.alpha).*exact.vPrev;
    exact.pa = pmc.alpha.*exact.pn + (1-pmc.alpha).*exact.pPrev;
    exact.DPDTa = pmc.alpha.*exact.DPDTn + (1-pmc.alpha).*exact.DPDTPrev;

    %analytical solution of pressure variation
    exact.k = exact.omega/acoustic.C;
    exact.phi = atan(1/(exact.k*0.5));
    exact.dp = fluid.dens.*acoustic.C.*exact.u.*(0.5./exact.r).*cos(exact.phi).*sin(time.*exact.omega - exact.k.*(exact.r - 0.5) + exact.phi);

    %Copying exact solution ( n + alpha ) to incompressible flow solution
    Sol.u(:,1) = exact.ua;
    Sol.u(:,2) = exact.va;
    Sol.p = exact.pa;
    Sol.q = -exact.DPDTa; % Source term

    % Nonlinear iterations start here
    for nLiter = 1:solver.nLIterMax

%         [Sol, NSnormIndicator, Force, Traction, Moment] = fluidDomain(solver, fluid, pmc, Sol, cnnF,... 
%                                              crdNew, elemType, ndof, nen, nElemF, BCFluid,...
%                                              fluidBCNodesU, fluidBCNodesV, ...
%                                              fluidBCNodesP,...
%                                              oisdNodes, r_0);

        [Sol, lpceNormIndicator] = lpceDomain(Sol, solver, acoustic, pmc, cnn, crd, elemType,...
                            nen, ndof, BCFluid, fluid,nElemF) ;

        if (lpceNormIndicator < solver.nLTol) %&& NSnormIndicator < solver.nLTol && norm(rk) < 1e-6)
            break;
        end
    end     
 
    fprintf('\n');
    
    % Copy current variables to previous variables
    Sol.uPrev = Sol.u ;
    Sol.uDotPrev = Sol.uDot ;
    Sol.pPrev = Sol.p ;
    Sol.aleDispPrev = Sol.aleDisp ;
    Sol.aleDispPrevPrev = Sol.aleDispPrev ;
    Sol.aleVelPrev = Sol.aleVel ;

    % Copy current variables to previous variables Acoustics 
    Sol.uaPrev = Sol.ua;
    Sol.vaPrev = Sol.va; 
    Sol.paPrev = Sol.pa; 
    Sol.rhoPrev = Sol.rho; 

    % Copy current variables to previous variables exact
    exact.vrPrev = exact.vrn;
    exact.uPrev = exact.un;
    exact.vPrev = exact.vn;
    exact.pPrev = exact.pn;
    exact.DPDTPrev = exact.DPDTn;

    %hydrodynamic density correction 
    exact.rho1 = fluid.dens.*exact.omega.*acoustic.M.*0.5^2.*cos(exact.omega*time)./(acoustic.C.*exact.r) +fluid.dens.*acoustic.M.^2.*0.5^4.*cos(2*exact.omega*time)./(4.*exact.r.^4);

    % plotting
    %X = (Sol.pa(fluid.radialNs) + exact.pn(fluid.radialNs) - (fluid.dens.*acoustic.C^2)/fluid.gamma) ./(fluid.dens.*acoustic.C^2);
    X = exact.rho1(fluid.radialNs) ./fluid.dens + Sol.pa(fluid.radialNs)./(fluid.dens.*acoustic.C^2);
    %X = (exact.rho1(fluid.radialNs) + Sol.rho(fluid.radialNs)) ./fluid.dens;
    Y = exact.dp(fluid.radialNs)./(fluid.dens.*acoustic.C^2);
    [pl(:,1),I] = sort(crd(fluid.radialNs,1)./0.5);
    pl(:,2) = X(I);
    pl(:,3) = Y(I);
    clear X Y

    figure(2);
    plot(pl(:,1) , pl(:,2)  ,'red');
    hold on
    plot(pl(:,1), pl(:,3), "blue");
    hold off
    title('Radial variation of Pa* at time step ',timeStep)
    legend('Numerical','actual')
    xlim([1 26])
    ylim([-0.14 0.14])
    xlabel('r/a')
    ylabel('Pa*')
%    F(timeStep) = getframe(gcf); 
 
    % Post-process the results
%     if (mod(timeStep,solver.intgOutFreq)==0)
%        fprintf(fileId1,'====timeStep: %d====\n',timeStep);
%        fprintf(fileId1,'Traction (total):\n');
%        fprintf(fileId1,'%e %e\n',Traction(1),Traction(2));
%        fprintf(fileId1,'Moment:\n');
%        fprintf(fileId1,'%e\n',Moment);
%     end
%      if (mod(timeStep,solver.intgOutFreq)==0)
%        fprintf(fileId2,'====timeStep: %d====\n',timeStep);
%        fprintf(fileId2,'Probe:\n');
%        fprintf(fileId2,'%d\n',probes);
%        fprintf(fileId2,'Displacement:\n');
%        fprintf(fileId2,'%e %e\n',Sol.aleDisp(probes,1),Sol.aleDisp(probes,2));
%     end
    
    if (timeStep >= solver.outFreq(1))
        outFreq = solver.outFreq(2) ;
    else
        outFreq = 1e10 ;
    end
    if (mod(timeStep,outFreq)==0)
        data = [crd(:,1)+Sol.aleDisp(:,1) crd(:,2)+Sol.aleDisp(:,2) crd(:,3) Sol.u Sol.p Sol.aleDisp ];
        filename = sprintf('%s/plt_files/%s.%d.plt',wrkDir,problemString,timeStep);
        fileId = fopen(filename,'w');

        FileTitle = 'Tecplot data';

        fprintf(fileId,' TITLE = \"%s\"\n',FileTitle);
        fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"Z\" \"U\", \"V\", \"p\", \"dispX\", \"dispY\", \"pDash\"\n');
    
        if (strcmp(elemType,'3Tri'))
            fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',(timeStep-1)*solver.dt,size(crd,1),size(cnn,1));
        elseif (strcmp(elemType,'4Quad'))
            fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',(timeStep-1)*solver.dt,size(crd,1),size(cnn,1));
        end
    
        fprintf(fileId,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n',data');

        if (strcmp(elemType,'3Tri'))
            fprintf(fileId,'%d %d %d\n',cnn');
        elseif (strcmp(elemType,'4Quad'))
            fprintf(fileId,'%d %d %d %d\n',cnn');
        end
        fclose(fileId);
        
        % Write vtk file
        datavtk = [crd(:,1)+Sol.aleDisp(:,1) crd(:,2)+Sol.aleDisp(:,2) crd(:,3) Sol.u Sol.p Sol.aleDisp Sol.pa./(fluid.dens.*acoustic.C^2) Sol.ua Sol.va Sol.rho];
        filenamevtk = sprintf('%s/vtk_files/%s.%d.vtk',wrkDir,problemString,timeStep);
        fileIdvtk = fopen(filenamevtk,'w');

        fprintf(fileIdvtk,'# vtk DataFile Version 2.0\n');
        fprintf(fileIdvtk,'TITLE = \"Vtk output\" \n');
        fprintf(fileIdvtk,'ASCII\n');
        fprintf(fileIdvtk,'DATASET UNSTRUCTURED_GRID\n');
        fprintf(fileIdvtk,'POINTS %d FLOAT\n',size(crd,1));
        fprintf(fileIdvtk,'%12.10f %12.10f %12.10f\n',[crd(:,1)+Sol.aleDisp(:,1) crd(:,2)+Sol.aleDisp(:,2) crd(:,3)]');
        if (strcmp(elemType,'3Tri'))
            nen = 3 ;
        elseif (strcmp(elemType,'4Quad'))
            nen = 4 ;
        end
        
        fprintf(fileIdvtk,'CELLS %d %d \n',size(cnn,1),size(cnn,1)*(nen+1));
        if (nen == 3)
            fprintf(fileIdvtk,'%d %d %d %d \n',[nen*ones(size(cnn,1),1) cnn(:,1)-1 cnn(:,2)-1 cnn(:,3)-1]');
        elseif (nen == 4)
            fprintf(fileIdvtk,'%d %d %d %d %d \n',[nen*ones(size(cnn,1),1) cnn(:,1)-1 cnn(:,2)-1 cnn(:,3)-1 cnn(:,4)-1]');
        end
        
        fprintf(fileIdvtk,'CELL_TYPES %d \n',size(cnn,1));
        if (nen == 3)
           fprintf(fileIdvtk,'%d \n',5*ones(size(cnn,1),1)); 
        elseif (nen == 4)
           fprintf(fileIdvtk,'%d \n',9*ones(size(cnn,1),1)); 
        end
        
        fprintf(fileIdvtk,'POINT_DATA %d \n',size(crd,1));
        
        fprintf(fileIdvtk,'SCALARS U FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.u(:,1));
        fprintf(fileIdvtk,'SCALARS V FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.u(:,2));
        fprintf(fileIdvtk,'SCALARS P FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.p);
        fprintf(fileIdvtk,'SCALARS dispX FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.aleDisp(:,1));
        fprintf(fileIdvtk,'SCALARS dispY FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.aleDisp(:,2));
        fprintf(fileIdvtk,'SCALARS Pa* FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.pa./(fluid.dens.*acoustic.C^2));
        fprintf(fileIdvtk,'SCALARS Ua FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.ua);
        fprintf(fileIdvtk,'SCALARS Va FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.va);
        fprintf(fileIdvtk,'SCALARS Rho FLOAT \n');
        fprintf(fileIdvtk,'LOOKUP_TABLE default \n');
        fprintf(fileIdvtk,'%12.10f \n',Sol.rho);
        fclose(fileIdvtk);
    end
    
end
fclose(fileId1);
fclose(fileId2);

%Storing the video
% video = VideoWriter("test.avi" , "Uncompressed AVI");
% video.FrameRate = 2;
% open(video)
% writeVideo(video,F);
% close(video)
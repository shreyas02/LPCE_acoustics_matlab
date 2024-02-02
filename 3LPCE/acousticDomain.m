function [Sol, AcNormIndicator] = acousticDomain(Sol, solver, acoustic, pmc, cnnF, crd, elemType,...
                            nenF, ndof, nElemF, BCFluid, acousticBCNodes)

% Solve acoustic equation
acousticBCValues = [acoustic.DirichletVal];

[Sol, AcNormIndicator] = acousticWaveEqn(Sol, solver, acoustic, pmc, cnnF, crd, elemType,...
                            nenF, ndof, nElemF, BCFluid, acousticBCNodes, ...
                            acousticBCValues);

end
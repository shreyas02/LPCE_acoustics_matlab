function [Sol, lpceNormIndicator] = lpceDomain(Sol, solver, acoustic, pmc, cnn, crd, elemType,...
                            nen, ndof, BCFluid , fluid,nElemF)

% Solve lpce equation

[Sol, lpceNormIndicator] = lpce(Sol, solver, acoustic, pmc, cnn, crd, elemType,...
                            nen, ndof, BCFluid , fluid,nElemF);
end
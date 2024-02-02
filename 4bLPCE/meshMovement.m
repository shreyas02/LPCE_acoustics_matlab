function [Sol, crdNew] = meshMovement(Sol, H, solver, crd, crdNew, cnn, elemType, nen, sourcePts, targetPts, BCInner, BCBound)

% Displace the nodes using RBF mapping
Sol.aleDisp(sourcePts,1) = [Sol.dispS(BCInner,1); zeros(size(BCBound,1),1)] ;
Sol.aleDisp(sourcePts,2) = [Sol.dispS(BCInner,2); zeros(size(BCBound,1),1)] ;
        
Sol.aleDisp(targetPts,1) = H*[Sol.dispS(BCInner,1); zeros(size(BCBound,1),1)];
Sol.aleDisp(targetPts,2) = H*[Sol.dispS(BCInner,2); zeros(size(BCBound,1),1)];
Sol.aleVel = (Sol.aleDisp - Sol.aleDispPrev)./solver.dt  ;
%Sol.aleVel = (3*Sol.aleDisp - 4*Sol.aleDispPrev + Sol.aleDispPrevPrev)./(2*solver.dt)  ;

%Sol.aleVel(sourcePts,1) = [Sol.uS(BCInner,1); zeros(size(BCBound,1),1)] ;
%Sol.aleVel(sourcePts,2) = [Sol.uS(BCInner,2); zeros(size(BCBound,1),1)] ;
%Sol.aleVelPrev(sourcePts,1) = [Sol.uSPrev(BCInner,1); zeros(size(BCBound,1),1)] ;
%Sol.aleVelPrev(sourcePts,2) = [Sol.uSPrev(BCInner,2); zeros(size(BCBound,1),1)] ;
                
clear crdNew ;
crdNew(:,1) = crd(:,1) + Sol.aleDisp(:,1);
crdNew(:,2) = crd(:,2) + Sol.aleDisp(:,2);
crdNew(:,3) = crd(:,3);
    
% Check mesh for negative Jacobians
checkMesh(crdNew, cnn, elemType, nen) ;

end
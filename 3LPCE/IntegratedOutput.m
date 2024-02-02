%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Calculate the integrated output values of the forces on the structure  %
%                                                                         %
%  This function evaluated the integral of the stress term on the surfa-  %
%  ce of the structure, i.e.,                                             %
%                                                                         %
%                   \int_\Gamma (\sigma . n) d\Gamma,                     %
%                                                                         %
%  where \sigma is the fluid stress, \sigma = -pI + \mu(grad U + grad U^T)%
%  p being the pressure, I is the Identity matrix, \mu is the fluid visc- %
%  osity and grad U is the gradient of velocity. n is the normal to the   %
%  structural surface.                                                    %
%                                                                         %
%  The integral is ealuated by reduced quadrature integration on the two- %
%  dimensional element just beside the surface.                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Traction, Moment] = IntegratedOutput(crd, Force, r_0)

% Tractions in X and Y directions
Traction(1) = sum(Force(:,1)) ;
Traction(2) = sum(Force(:,2)) ;

% Z-Moment
Moment = Force(:,2).*(crd(:,1) - r_0(1)) - Force(:,1).*(crd(:,2) - r_0(2)) ;
Moment = sum(Moment(:)) ;

end

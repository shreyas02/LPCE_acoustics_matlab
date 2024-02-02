function exact = incomp(crd , time);

%Mach number 
mach = 0.08;
cs = 330;

%fluid properties
rho = 1;
g = 1.4; %gamma of fluid 

%coordinates
x = crd(:,1);
y = crd(:,2);
r = sqrt(x.^2+y.^2);
theta = atan(y./x);
z = x + 1i.*y;
r0 = 1; %(assumed)

%timesteps 

%calculation of gamma and omega
gamma = mach*4*pi*r0*cs;
omega = gamma/(4*pi*r0^2);

b = r0*exp(1i*omega*time);
%calculation of the velocity feild 
test = (gamma/(pi*1i)).*(z./(z.^2 + b^2));
u = real(test);
v = -imag(test);
clear test
%calculation of the velocity magnitude
X = real(z.^2 - b^2);
Y = imag(z.^2 - b^2);
R = sqrt(X.^2 +Y.^2);
Q = (gamma.*r)./(pi.*R);
%calculation of phiT
phiT = ((-omega.*gamma)./pi).*real((b^2)./(z.^2 - b^2));
%calculation of p0
p0 = (rho * cs * cs)/g;
%incompressible pressure
p = p0 - rho.*phiT - (rho.*Q.^2)./2;
%velocity time derivatives
test = ((2*omega*gamma)/(pi)).*((z.*b.^2)./((z.^2-b.^2).^2));
ut = real(test);
vt = -imag(test);
clear test
%velocity x derivative
test = ((-gamma)/(pi.*1i)).*((z.^2 + b.^2)./((z.^2-b.^2).^2));
ux = real(test);
vx = -imag(test);
clear test
%velocity y derivative
uy = vx;
vy = -ux;
%calculation of phiTT
phiTT = ((2*omega^2*gamma)/(pi))*imag((z.^2.*b.^2)./((z.^2-b.^2).^2));
%p t derivative
pt = -rho.*(phiTT + u.*ut + v.*vt);
%p x derivative
px = -rho.*(ut + u.*ux + v.*vx);
%p y derivative
py = -rho.*(vt + u.*uy + v.*vy);
%total p derivative
DPDT = pt + u.*px + v.*py;
%calculation of acoustic pressure
k = (2*omega)/cs;
H = besselh(2,2,k.*r);
pa = ((1i*rho*gamma^4)/(64*pi^3*r0^4*cs^2)).*H.*exp(2i.*(omega*time - theta)) ;
pa = real(pa);
%copying variables to exact variable

exact.v = v;
exact.u = u;
exact.p = p;
exact.dpdt = DPDT;
exact.pa = pa;
exact.r = r;






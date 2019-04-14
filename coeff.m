* Matlab code for Cp calculation  
input: unscaled moment of one blade

output: power coefficient of a 3-blades wind/tidal turbine 

function cp=coeff(m)
% v : velocity
% m : unscaled moment of one blade
% rho : density
% cp : power coefficient of 3 blades
% omega: angular velocity, rad/s
v = 0.6;
A = 0.166;
omega = 11.087;
rho = 998.2;
cp=(3*m*omega)/(0.5*rho*v*v*v*A)

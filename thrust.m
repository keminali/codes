function ct=thrust(t)
% v : velocity
% t: unscaled thrust of one blade
% rho : density
% ct : thrust coefficient of 3 blades
% omega: angular velocity, rad/s
v = 0.6;
A = 0.166;
rho = 998.2;
ct=(3*t)/(0.5*rho*v*v*A)

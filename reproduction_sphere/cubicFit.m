%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cubic-fit interpolation
% Interpolate the states and controlled variable
% Input: 
%   s0: Initial state - theta0, phi0, d(theta)/dt(0), d(phi)/dt(0)
%   sf: Final state - thetaf, phif, d(theta)/dt(f), d(phi)/dt(f)
%   t0: Initial time
%   tf: Final time
%   M1,M2: C-W matrix
%   rho: radial length
% Output:
%   J: Optimizing index
%   u: controlled variable
%   state: theta, phi, d(theta)/dt, d(phi)/dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J,u,SphericalState4D] = cubicFit(SphericalState4D0, SphericalState4Df, t0, tf, M1, M2, rho)
% F Matrix
F = [1, t0, t0^2, t0^3;
     0, 1,  2*t0, 3*t0^2;
     1, tf, tf^2, tf^3;
     0, 1,  2*tf, 3*tf^2];

% Initial state
theta0 = SphericalState4D0(1);
phi0 = SphericalState4D0(2);
vTheta0 = SphericalState4D0(3);
vPhi0 = SphericalState4D0(4);

% Final state
thetaf = SphericalState4Df(1);
phif = SphericalState4Df(2);
vThetaf = SphericalState4Df(3);
vPhif = SphericalState4Df(4);

% Interpolation coefficient
bTheta = [theta0, vTheta0, thetaf, vThetaf]';
bPhi = [phi0, vPhi0, phif, vPhif]';
c = F \ bTheta;
d = F \ bPhi;

% Interpolation and get the optimizing index
syms t
assume(t >= 0)

% Interpolate the state and acceleration
% State
SphericalState4D(1) = c(4) * t^3 + c(3) * t^2 + c(2) * t + c(1);
SphericalState4D(2) = d(4) * t^3 + d(3) * t^2 + d(2) * t + d(1);
SphericalState4D(3) = 3 * c(4) * t^2 + 2 * c(3) * t + c(2);
SphericalState4D(4) = 3 * d(4) * t^2 + 2 * d(3) * t + d(2);


% Acceleration
SphericalAcc(1) = 6 * c(4) * t + 2 * c(3);
SphericalAcc(2) = 6 * d(4) * t + 2 * d(3);

% Variable transform to make the code cleaner
vTheta = SphericalState4D(3);
vPhi = SphericalState4D(4);
aTheta = SphericalAcc(1);
aPhi = SphericalAcc(2);
sTheta = sin(SphericalState4D(1));
cTheta = cos(SphericalState4D(1));
sPhi = sin(SphericalState4D(2));
cPhi = cos(SphericalState4D(2));

% Cartesian position and velocity
p = rho .* [cTheta .* sPhi;
            sTheta .* sPhi;
            cPhi];                          % Position
v = rho .* [-vTheta .* sTheta .* sPhi + vPhi .* cTheta .* cPhi;             % Velocity
            vTheta .* cTheta .* sPhi + vPhi .* sTheta .* cPhi;   
            -vPhi .* sPhi];

a = rho .* [-aTheta .* sTheta .* sPhi - vTheta.^2 .* cTheta .* sPhi ...     % Acceleration
           - 2 .* vTheta .* vPhi .* sTheta .* cPhi + aPhi .* cTheta .* cPhi ...
           - vTheta.^2 .* cTheta .* sPhi,...
            aTheta .* cTheta .* sPhi - vTheta.^2 .* sTheta .* sPhi ...
           + 2 .* vTheta .* vPhi .* cTheta .* cPhi + aPhi .* sTheta .* cPhi ...
           - vPhi.^2 .* sTheta .* sPhi, ...
           -aPhi .* sPhi - vPhi.^2 .* cPhi]';

% Controlled variable
u = a - M1 * p - M2 * v;

% The optimal index
E = 0.5 * (u' * u);
energy = @(x)double(subs(E, t, x));
J = gaussLegendre5_comp(energy, t0, tf, 10);
end


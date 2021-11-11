clc
clear all
close all

%% Needed scripts path
addpath(strcat(pwd, '\timeConversion'))

%% DEFINE SETTINGS STRUCTURE
% Create a Simulink object to use settings struct in Simulink
settings = Simulink.Parameter;
settings.Value = struct;
settings.CoderInfo.StorageClass = 'ExportedGlobal';

%% SETUP
%%% Simulation configuration
settings.Value.Nperiod = 1;             % [-] Number of nominal period
settings.Value.quat = true;             % [-] True if kinematics are computed with quaternions

%%% Orbit perturbation models, followed by model constants
settings.Value.orbitPert.J2 = true;     % [-] True if J2 orbit perturbation is computed
settings.Value.orbitPert.J2value = astroConstants(9);

%%% Starting orbit parameters (INERTIAL FRAME)
settings.Value.mu = astroConstants(13); % [km^3/s^2] Earth's gravitational parameter
settings.Value.Re = astroConstants(23); % [km] Earth's mean radius
settings.Value.r0 = [26578.137; 0; 0];  % [km] Initial position in reference frame
settings.Value.v0 = [0; 2.221; 1];      % [km/s] Initial velocity in reference frame
% Create initial state vector
settings.Value.Y0 = [settings.Value.r0; settings.Value.v0];

%%% Inertias
settings.Value.Ix = 0.01;               % [kg m^2] Inertia moment along X axis
settings.Value.Iy = 0.05;               % [kg m^2] Inertia moment along Y axis
settings.Value.Iz = 0.07;               % [kg m^2] Inertia moment along Z axis
% Create the inertia matrix and its inverse matrix
settings.Value.J = diag([settings.Value.Ix settings.Value.Iy settings.Value.Iz]);
settings.Value.invJ = inv(settings.Value.J);

%%% Initial angular velocities (BODY FRAME)
settings.Value.wx0 = 0.001;              % [rad/s] Initial X axis angular velocity
settings.Value.wy0 = 0.005;              % [rad/s] Initial Y axis angular velocity
settings.Value.wz0 = 0.001;              % [rad/s] Initial Z axis angular velocity
% Create the intial angular velocity vector
settings.Value.omega0 = [settings.Value.wx0 settings.Value.wy0 settings.Value.wz0]';

%%% Initial direction cosines matrix
settings.Value.DCM0 = eye(3);

%%% Initial quaternion
settings.Value.q0 = [1 0 0 0]';

%% PRE SIMULATION CALCULATIONS
r0 = settings.Value.r0;
v0 = settings.Value.v0;
mu = settings.Value.mu;
a = 1/( 2/norm(r0) - dot(v0,v0)/mu );    % [km] Semi-major axis
Tperiod = 2*pi*sqrt( a^3/mu );           % [s] Orbital period

% Get Julian day 2000
t = datetime('now');
[y, M, d] = ymd(t);
[h, m, s] = hms(t);
date = [y, M, d, h, m, s];
time_mjd2000 = date2mjd2000(date);

% Get Keplerian parameters of the Earth
[kep_Earth, ksun] = uplanet(time_mjd2000, 3);

%% Get Cartesian position and velocity of the Earth
% Planetary orbital elements are restituited in a Sun-centred ecliptic system.
[rrE, vvE] = par2car(kep_Earth, ksun);
settings.Value.EarthY0 = [rrE; vvE];

%% Get Cartesian position and velocity of the Moon
% Position vector of the Moon in cartesian coordinates, expressed in the
% Geocentric Equatorial Reference Frame (IAU-76/FK5 J2000, mean equator,
%   mean equinox frame)
[rrM, vvM] = ephMoon(time_mjd2000);
settings.Value.MoonY0 = [rrM; vvM];


%% RUN THE SIMULATION
% Set the simulation time
settings.Value.Tsim = settings.Value.Nperiod * Tperiod;

% Simulation run
simOut = sim('model.slx', 'SrcWorkspace', 'current');

% Retrieve output data
T = simOut.tout;            % [s] Simulation time
R = simOut.Y.Data;
DCM_I2L = simOut.DCM_I2L.Data;
DCM_I2B = simOut.DCM_I2B.Data;
DCM_L2B = simOut.DCM_L2B.Data;
quat = simOut.quat.Data;

%% PLOT
figure;
quiver3(0,0,0,10000,0,0); grid on; hold on
axis equal
quiver3(0,0,0,0,10000,0);
quiver3(0,0,0,0,0,10000);

% Local
xl_vec = 5000 .* DCM_I2L(:,:,1)'*[1 0 0]';
xl = quiver3(R(1,1),R(1,2),R(1,3),xl_vec(1),xl_vec(2),xl_vec(3));
yl_vec = 5000 .* DCM_I2L(:,:,1)'*[0 1 0]';
yl = quiver3(R(1,1),R(1,2),R(1,3),yl_vec(1),yl_vec(2),yl_vec(3));
zl_vec = 5000 .* DCM_I2L(:,:,1)'*[0 0 1]';
zl = quiver3(R(1,1),R(1,2),R(1,3),zl_vec(1),zl_vec(2),zl_vec(3));

% Body
xb_vec = 5000 .* DCM_I2B(:,:,1)'*[1 0 0]';
xb = quiver3(R(1,1),R(1,2),R(1,3),xb_vec(1),xb_vec(2),xb_vec(3));
yb_vec = 5000 .* DCM_I2B(:,:,1)'*[0 1 0]';
yb = quiver3(R(1,1),R(1,2),R(1,3),yb_vec(1),yb_vec(2),yb_vec(3));
zb_vec = 5000 .* DCM_I2B(:,:,1)'*[0 0 1]';
zb = quiver3(R(1,1),R(1,2),R(1,3),zb_vec(1),zb_vec(2),zb_vec(3));

plot3(R(:,1), R(:,2), R(:,3), 'r', 'LineWidth', 1)

for i = 2:100:length(T)
    % Local
    delete(xl)
    delete(yl)
    delete(zl)
    xl_vec = 5000 .* DCM_I2L(:,:,i)'*[1 0 0]';
    xl = quiver3(R(i,1),R(i,2),R(i,3),xl_vec(1),xl_vec(2),xl_vec(3),'k');
    yl_vec = 5000 .* DCM_I2L(:,:,i)'*[0 1 0]';
    yl = quiver3(R(i,1),R(i,2),R(i,3),yl_vec(1),yl_vec(2),yl_vec(3),'k');
    zl_vec = 5000 .* DCM_I2L(:,:,i)'*[0 0 1]';
    zl = quiver3(R(i,1),R(i,2),R(i,3),zl_vec(1),zl_vec(2),zl_vec(3),'k');

    % Body
    delete(xb)
    delete(yb)
    delete(zb)
    xb_vec = 5000 .* DCM_I2B(:,:,i)'*[1 0 0]';
    xb = quiver3(R(i,1),R(i,2),R(i,3),xb_vec(1),xb_vec(2),xb_vec(3),'b');
    yb_vec = 5000 .* DCM_I2B(:,:,i)'*[0 1 0]';
    yb = quiver3(R(i,1),R(i,2),R(i,3),yb_vec(1),yb_vec(2),yb_vec(3),'r');
    zb_vec = 5000 .* DCM_I2B(:,:,i)'*[0 0 1]';
    zb = quiver3(R(i,1),R(i,2),R(i,3),zb_vec(1),zb_vec(2),zb_vec(3),'y');
    drawnow limitrate
end


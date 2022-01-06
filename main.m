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

%% *************************** USER HERE **********************************
%%% SIMULATION CONFIGURATION
settings.Value.date0   = [2021 12 20 12 0 0];
settings.Value.Nperiod = 1;             % [-] Number of nominal period

%%% PERTURBATION MODELS
settings.Value.orbitPert.J2 = true;     % [-] True if J2 orbit perturbation is computed

%%% Disturbances
settings.Value.r_CG = [0.02 0 0];
settings.Value.CD = 2.2;
settings.Value.rho = 1.454*1e-13;
settings.Value.dec = 10;

%%% STARTING ORBIT PARAMETERS (INERTIAL FRAME)
a  = 6971;                             % [km] Orbit semi-major axis
e  = 0.01;                            % [-] Orbit eccentricity
i  = deg2rad(10);                 % [rad] Orbit innclination
OM = deg2rad(80);                  % [rad] Orbit RAAN
om = deg2rad(123);                  % [rad] Orbit pericenter anomaly
th = deg2rad(0);                                 % [rad] Orbit true anomaly

%%% INERTIAS
settings.Value.Ix = 1.009;               % [kg m^2] Inertia moment along X axis
settings.Value.Iy = 0.251;               % [kg m^2] Inertia moment along Y axis
settings.Value.Iz = 0.961;               % [kg m^2] Inertia moment along Z axis

%%% INITIAL ANGULAR VELOCITIES (BODY FRAME)
settings.Value.wx0 = deg2rad(-12);             % [rad/s] Initial X axis angular velocity
settings.Value.wy0 = deg2rad(-8);             % [rad/s] Initial Y axis angular velocity
settings.Value.wz0 = deg2rad(7);             % [rad/s] Initial Z axis angular velocity

%%% Initial direction cosines matrix
settings.Value.DCM0 = [1 0 0; 0 1 0; 0 0 1];

%%% Initial quaternion
settings.Value.q0 = dcm2quat([1 0 0; 0 1 0; 0 0 1]);

%%% SENSORS
settings.Value.SunSensorAccuracy = 0.3;
settings.Value.SunSensorSampleRate = 50;
settings.Value.EarthSensorAccuracy = 0.25;
settings.Value.EarthSensorSampleRate = 100;
settings.Value.GyroscopeARW = 0.2;
settings.Value.GyroscopeRRW = 0.3;
settings.Value.GyroscopeSampleRate = 1000;

settings.Value.MagnetoTorquerMaxDipole = 1.9;
settings.Value.MagnetoTorquerDipoleError = 0.005;
settings.Value.MagnetoTorquerSampleRate = 10;

settings.Value.ObserverLw = 0.05;
settings.Value.ObserverLd = 5.8e-4;
settings.Value.ObserverMd0 = [0 0 0]';

settings.Value.detumblingOnlyMagnMaxTime = 1500;
settings.Value.detumblingDetumblingMaxTime = 2500;
settings.Value.slewManeuverMaxTime = 3500;
settings.Value.detumblingkb = 1e7;
settings.Value.detumblingkp = 0.05;
settings.Value.slewMank1 = 0.1;
settings.Value.slewMank2 = 0.001;
settings.Value.earthPointk1 = 0.5;
settings.Value.earthPointk2 = 0.05;

%%% Reaction Wheels
settings.Value.IRW = 1e-5;              % [kg*m^2] Reaction Wheels' Moment of Inertia
settings.Value.MaxRPM = 8000;           % Reaction Wheels' Maximum RPM Allowed
settings.Value.MaxRPM_Rate = 100;       % Reaction Wheels' Maximum RPM Rate Allowed
settings.Value.RPMAccuracy = 5;         % Reaction Wheels' RPM Accuracy
settings.Value.RW_sr = 1;               % Reaction Wheels' Sampling Rate

settings.Value.alpha1StateObserver = -0.6;
settings.Value.alpha1StateObserverGyro = -0.4;


%% ***** FROM NOW ON DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING *****
%%% PERTURBATION MODELS
settings.Value.orbitPert.J2value = astroConstants(9);

%%% STARTING ORBIT PARAMETERS (INERTIAL FRAME)
settings.Value.muE = astroConstants(13); % [km^3/s^2] Earth's gravitational parameter
settings.Value.muS = astroConstants(4);  % [km^3/s^2] Sun's gravitational parameter
settings.Value.Re = astroConstants(23); % [km] Earth's mean radius
% Calculates Initial orbit position and velocity
[rr0, vv0] = par2car([a, e, i, OM, om, th], settings.Value.muE);
settings.Value.r0 = rr0;                % [km] Initial position in reference frame
settings.Value.v0 = vv0;                % [km/s] Initial velocity in reference frame

% Create initial state vector
settings.Value.Y0 = [settings.Value.r0; settings.Value.v0];
kepEarth = uplanet(date2mjd2000(settings.Value.date0), 3);
[rrE, vvE] = par2car(kepEarth, settings.Value.muS);
settings.Value.Y0E = [rrE; vvE];

%%% INERTIAS
% Create the inertia matrix and its inverse matrix
settings.Value.J = diag([settings.Value.Ix settings.Value.Iy settings.Value.Iz]);
settings.Value.invJ = inv(settings.Value.J);

%%% INITIAL AGNGULAR VELOCITIES
% Create the intial angular velocity vector
settings.Value.omega0 = [settings.Value.wx0 settings.Value.wy0 settings.Value.wz0]';



%% RUN THE SIMULATION
% Set the simulation time
Tperiod = 2*pi * sqrt(a^3/settings.Value.muE);
% settings.Value.Tsim = settings.Value.Nperiod * Tperiod/5;
settings.Value.Tsim = Tperiod;

%%
% Simulation run
simOut = sim('model.slx', 'SrcWorkspace', 'current');

%%
Retrieve output data
T = simOut.tout;            % [s] Simulation time
R = simOut.Y.Data;
DCM_I2L = simOut.DCM_I2L.Data;
DCM_I2B = simOut.DCM_I2B.Data;
DCM_L2B = simOut.DCM_L2B.Data;
quaternions = simOut.quaternions.Data;
omegaReal = simOut.omegaReal.Data;
disturb = simOut.Disturbances.Data;
eclipse = simOut.eclipse.Data;
attErr = simOut.attitudeEstimate.Data;
estMD = simOut.estMD.Data;
magnetoM = simOut.MagnetoM.Data;
RWsM = simOut.RWsM.Data;
Mtot = simOut.Mtot.Data;

%% detumbling only magnetotorquers
index = find(T >= 0 & T < 1500);
plot(index,quaternions);
plot(index,omegaReal);
plot(index,disturb);
plot(index,eclipse);
plot(index,attErr);
plot(index,estMD);
plot(index,magnetoM);
plot(index,RWsM);
plot(index,Mtot);

%% second phase detumbling
index = find(T >1500 & T < 2500);
plot(index,quaternions);
plot(index,omegaReal);
plot(index,disturb);
plot(index,eclipse);
plot(index,attErr);
plot(index,estMD);
plot(index,magnetoM);
plot(index,RWsM);
plot(index,Mtot);

%% slew maneuver
index = find(T > 2500 & T < 3500);
plot(index,quaternions);
plot(index,omegaReal);
plot(index,disturb);
plot(index,eclipse);
plot(index,attErr);
plot(index,estMD);
plot(index,magnetoM);
plot(index,RWsM);
plot(index,Mtot);

%% earth poiting maneuver
index = find(T >3500 & T < Tperiod);
plot(index,quaternions);
plot(index,omegaReal);
plot(index,disturb);
plot(index,eclipse);
plot(index,attErr);
plot(index,estMD);
plot(index,magnetoM);
plot(index,RWsM);
plot(index,Mtot);

%% PLOT
% figure;
% quiver3(0,0,0,10000,0,0); grid on; hold on
% axis equal
% quiver3(0,0,0,0,10000,0);
% quiver3(0,0,0,0,0,10000);
% 
% % Local
% xl_vec = 5000 .* DCM_I2L(:,:,1)'*[1 0 0]';
% xl = quiver3(R(1,1),R(1,2),R(1,3),xl_vec(1),xl_vec(2),xl_vec(3), 'k');
% yl_vec = 5000 .* DCM_I2L(:,:,1)'*[0 1 0]';
% yl = quiver3(R(1,1),R(1,2),R(1,3),yl_vec(1),yl_vec(2),yl_vec(3), 'k');
% zl_vec = 5000 .* DCM_I2L(:,:,1)'*[0 0 1]';
% zl = quiver3(R(1,1),R(1,2),R(1,3),zl_vec(1),zl_vec(2),zl_vec(3), 'k');
% 
% % Body
% xb_vec = 5000 .* DCM_I2B(:,:,1)'*[1 0 0]';
% xb = quiver3(R(1,1),R(1,2),R(1,3),xb_vec(1),xb_vec(2),xb_vec(3), 'b');
% yb_vec = 5000 .* DCM_I2B(:,:,1)'*[0 1 0]';
% yb = quiver3(R(1,1),R(1,2),R(1,3),yb_vec(1),yb_vec(2),yb_vec(3), 'r');
% zb_vec = 5000 .* DCM_I2B(:,:,1)'*[0 0 1]';
% zb = quiver3(R(1,1),R(1,2),R(1,3),zb_vec(1),zb_vec(2),zb_vec(3), 'y');
% 
% plot3(R(:,1), R(:,2), R(:,3), '--r', 'LineWidth', 1)
% 
% for i = 2:300:length(T)
%     % Local
%     delete(xl)
%     delete(yl)
%     delete(zl)
%     xl_vec = 5000 .* DCM_I2L(:,:,i)'*[1 0 0]';
%     xl = quiver3(R(i,1),R(i,2),R(i,3),xl_vec(1),xl_vec(2),xl_vec(3),'k');
%     yl_vec = 5000 .* DCM_I2L(:,:,i)'*[0 1 0]';
%     yl = quiver3(R(i,1),R(i,2),R(i,3),yl_vec(1),yl_vec(2),yl_vec(3),'k');
%     zl_vec = 5000 .* DCM_I2L(:,:,i)'*[0 0 1]';
%     zl = quiver3(R(i,1),R(i,2),R(i,3),zl_vec(1),zl_vec(2),zl_vec(3),'k');
% 
%     % Body
%     delete(xb)
%     delete(yb)
%     delete(zb)
%     xb_vec = 5000 .* DCM_I2B(:,:,i)'*[1 0 0]';
%     xb = quiver3(R(i,1),R(i,2),R(i,3),xb_vec(1),xb_vec(2),xb_vec(3),'b');
%     yb_vec = 5000 .* DCM_I2B(:,:,i)'*[0 1 0]';
%     yb = quiver3(R(i,1),R(i,2),R(i,3),yb_vec(1),yb_vec(2),yb_vec(3),'r');
%     zb_vec = 5000 .* DCM_I2B(:,:,i)'*[0 0 1]';
%     zb = quiver3(R(i,1),R(i,2),R(i,3),zb_vec(1),zb_vec(2),zb_vec(3),'y');
%     drawnow limitrate
% end


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
a  = 6971;                                       % [km] Orbit semi-major axis
e  = 0.01;                                       % [-] Orbit eccentricity
i  = deg2rad(10);                                % [rad] Orbit innclination
OM = deg2rad(80);                                % [rad] Orbit RAAN
om = deg2rad(123);                               % [rad] Orbit pericenter anomaly
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
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Retrieve output data
T = simOut.tout;            % [s] Simulation time
% R = simOut.Y.Data;
DCM_I2L = simOut.DCM_I2L.Data;
DCM_I2B = simOut.DCM_I2B.Data;
DCM_L2B = simOut.DCM_L2B.Data;
quaternions = simOut.quaternions.Data;
omegaReal = simOut.omegaReal.Data;
omegaEstim = simOut.omega_estim.Data;
disturb = simOut.Disturbances.Data;
eclipse = simOut.eclipse.Data;
attErr = simOut.attitudeEstimate.Data;
estMD = simOut.estMD.Data;
magnetoM = simOut.MagnetoM.Data;
Mreal = simOut.Mtot.Data;
Gyro_end = simOut.gyro_end.Data;
Mc_RW = simOut.Mc_RW.Data;
Mid = simOut.Mc_magneto.Data;
pointing_err = simOut.pointing_err.Data;
RWsM = simOut.RWsM.Data;
Dmagn = simOut.Dmagn.Data;
a = 6971;
%% detumbling only magnetotorquers
indexDet = find(T >= 0 & T < settings.Value.detumblingDetumblingMaxTime);
indexOnlyMagn = find(T >= 0 & T < settings.Value.detumblingOnlyMagnMaxTime);
indexDetRW = find(T >= settings.Value.detumblingOnlyMagnMaxTime & T < settings.Value.detumblingDetumblingMaxTime);

%% DETUMBLING PLOTS
%%% omega Real only magnetotoruqers
figure
plot(T(indexOnlyMagn), omegaReal(indexOnlyMagn, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexOnlyMagn), omegaReal(indexOnlyMagn, 2), 'LineWidth', 1);
plot(T(indexOnlyMagn), omegaReal(indexOnlyMagn, 3), 'LineWidth', 1);
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
xlabel('t [s]'); ylabel('$\omega$ [$\frac{rad}{s}$]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% omega Real with RWs also
figure
plot(T(indexDetRW), omegaReal(indexDetRW, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexDetRW), omegaReal(indexDetRW, 2), 'LineWidth', 1);
plot(T(indexDetRW), omegaReal(indexDetRW, 3), 'LineWidth', 1);
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
xlabel('t [s]'); ylabel('$\omega$ [$\frac{rad}{s}$]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% quaternions
figure
plot(T(indexDet), quaternions(indexDet, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexDet), quaternions(indexDet, 2), 'LineWidth', 1);
plot(T(indexDet), quaternions(indexDet, 3), 'LineWidth', 1);
plot(T(indexDet), quaternions(indexDet, 4), 'LineWidth', 1);
legend('$q_0$', '$q_1$', '$q_2$', '$q_3$')
xlabel('t [s]'); ylabel('q [-]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% pointing error
figure
plot(T(indexDet), pointing_err(indexDet), 'LineWidth', 1); grid on; hold on;
xlabel('t [s]'); ylabel('Error [deg]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% Real control momentum
figure
plot(T(indexDet), Mreal(indexDet, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexDet), Mreal(indexDet, 2), 'LineWidth', 1);
plot(T(indexDet), Mreal(indexDet, 3), 'LineWidth', 1);
legend('$M_{c, x}$', '$M_{c, y}$', '$M_{c, z}$')
xlabel('t [s]'); ylabel('$M_c$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% Ideal - real control moment
figure
T1 = T(indexDetRW);
Mreal1 = Mreal(indexDetRW, :);
Mid1 = Mid(indexDetRW, :);
plot(T1(2:end), -Mreal1(2:end, 1)+Mid1(2:end, 1), 'LineWidth', 1); grid on; hold on;
plot(T1(2:end), -Mreal1(2:end, 2)+Mid1(2:end, 2), 'LineWidth', 1);
plot(T1(2:end), -Mreal1(2:end, 3)+Mid1(2:end, 3), 'LineWidth', 1);
legend('$\delta M_{c, x}$', '$\delta M_{c, y}$', '$\delta M_{c, z}$')
xlabel('t [s]'); ylabel('${M_{c}}^{id}$ - $M_c$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% RW torque
figure
plot(T(indexDet), RWsM(indexDet, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexDet), RWsM(indexDet, 2), 'LineWidth', 1);
plot(T(indexDet), RWsM(indexDet, 3), 'LineWidth', 1);
xlabel('t [s]'); ylabel('$M_r$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% magnetotorquer Dipole
figure
plot(T(indexDet), Dmagn(indexDet, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexDet), Dmagn(indexDet, 2), 'LineWidth', 1);
plot(T(indexDet), Dmagn(indexDet, 3), 'LineWidth', 1);
legend('$D_{x}$', '$D_{y}$', '$D_{z}$')
xlabel('t [s]'); ylabel('D [A$m^2$]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)


%% SLEW MANEUVER PLOTS
indexSlew = find(T >= settings.Value.detumblingDetumblingMaxTime & T < settings.Value.slewManeuverMaxTime);

%%% omega
figure
plot(T(indexSlew), omegaReal(indexSlew, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexSlew), omegaReal(indexSlew, 2), 'LineWidth', 1);
plot(T(indexSlew), omegaReal(indexSlew, 3), 'LineWidth', 1);
n = sqrt(398600/a^3);
plot(T(indexSlew), n*ones(length(T(indexSlew)), 1), '--k', 'LineWidth', 1)
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
xlabel('t [s]'); ylabel('$\omega$ [$\frac{rad}{s}$]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% quaternions
figure
plot(T(indexSlew), quaternions(indexSlew, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexSlew), quaternions(indexSlew, 2), 'LineWidth', 1);
plot(T(indexSlew), quaternions(indexSlew, 3), 'LineWidth', 1);
plot(T(indexSlew), quaternions(indexSlew, 4), 'LineWidth', 1);
legend('$q_0$', '$q_1$', '$q_2$', '$q_3$')
xlabel('t [s]'); ylabel('q [-]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% pointing error
figure
plot(T(indexSlew), pointing_err(indexSlew), 'LineWidth', 1); grid on; hold on;
xlabel('t [s]'); ylabel('Error [deg]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% Real control momentum
figure
plot(T(indexSlew), Mreal(indexSlew, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexSlew), Mreal(indexSlew, 2), 'LineWidth', 1);
plot(T(indexSlew), Mreal(indexSlew, 3), 'LineWidth', 1);
legend('$M_{c, x}$', '$M_{c, y}$', '$M_{c, z}$')
xlabel('t [s]'); ylabel('$M_c$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% Ideal - real control moment
figure
T1 = T(indexSlew);
Mreal1 = Mreal(indexSlew, :);
Mid1 = Mid(indexSlew, :);
plot(T1(2:end), -Mreal1(2:end, 1)+Mid1(2:end, 1), 'LineWidth', 1); grid on; hold on;
plot(T1(2:end), -Mreal1(2:end, 2)+Mid1(2:end, 2), 'LineWidth', 1);
plot(T1(2:end), -Mreal1(2:end, 3)+Mid1(2:end, 3), 'LineWidth', 1);
legend('$\delta M_{c, x}$', '$\delta M_{c, y}$', '$\delta M_{c, z}$')
xlabel('t [s]'); ylabel('${M_{c}}^{id}$ - $M_c$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% RW torque
figure
plot(T(indexSlew), RWsM(indexSlew, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexSlew), RWsM(indexSlew, 2), 'LineWidth', 1);
plot(T(indexSlew), RWsM(indexSlew, 3), 'LineWidth', 1);
xlabel('t [s]'); ylabel('$M_r$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% magnetotorquer Dipole
figure
plot(T(indexSlew), Dmagn(indexSlew, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexSlew), Dmagn(indexSlew, 2), 'LineWidth', 1);
plot(T(indexSlew), Dmagn(indexSlew, 3), 'LineWidth', 1);
legend('$D_{x}$', '$D_{y}$', '$D_{z}$')
xlabel('t [s]'); ylabel('D [A$m^2$]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)



%% EARTH POINT PLOTS
indexPoint = find(T >= settings.Value.slewManeuverMaxTime);

%%% omega
figure
plot(T(indexPoint), omegaReal(indexPoint, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexPoint), omegaReal(indexPoint, 2), 'LineWidth', 1);
plot(T(indexPoint), omegaReal(indexPoint, 3), 'LineWidth', 1);
n = sqrt(398600/a^3);
plot(T(indexPoint), n*ones(length(T(indexPoint)), 1), '--k', 'LineWidth', 1)
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
xlabel('t [s]'); ylabel('$\omega$ [$\frac{rad}{s}$]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% quaternions
figure
plot(T(indexPoint), quaternions(indexPoint, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexPoint), quaternions(indexPoint, 2), 'LineWidth', 1);
plot(T(indexPoint), quaternions(indexPoint, 3), 'LineWidth', 1);
plot(T(indexPoint), quaternions(indexPoint, 4), 'LineWidth', 1);
legend('$q_0$', '$q_1$', '$q_2$', '$q_3$')
xlabel('t [s]'); ylabel('q [-]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% pointing error
figure
plot(T(indexPoint), pointing_err(indexPoint), 'LineWidth', 1); grid on; hold on;
xlabel('t [s]'); ylabel('Error [deg]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% Real control momentum
figure
plot(T(indexPoint), Mreal(indexPoint, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexPoint), Mreal(indexPoint, 2), 'LineWidth', 1);
plot(T(indexPoint), Mreal(indexPoint, 3), 'LineWidth', 1);
legend('$M_{c, x}$', '$M_{c, y}$', '$M_{c, z}$')
xlabel('t [s]'); ylabel('$M_c$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% Ideal - real control moment
figure
T1 = T(indexPoint);
Mreal1 = Mreal(indexPoint, :);
Mid1 = Mid(indexPoint, :);
plot(T1(2:end), -Mreal1(2:end, 1)+Mid1(2:end, 1), 'LineWidth', 1); grid on; hold on;
plot(T1(2:end), -Mreal1(2:end, 2)+Mid1(2:end, 2), 'LineWidth', 1);
plot(T1(2:end), -Mreal1(2:end, 3)+Mid1(2:end, 3), 'LineWidth', 1);
legend('$\delta M_{c, x}$', '$\delta M_{c, y}$', '$\delta M_{c, z}$')
xlabel('t [s]'); ylabel('${M_{c}}^{id}$ - $M_c$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% RW torque
figure
plot(T(indexPoint), RWsM(indexPoint, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexPoint), RWsM(indexPoint, 2), 'LineWidth', 1);
plot(T(indexPoint), RWsM(indexPoint, 3), 'LineWidth', 1);
xlabel('t [s]'); ylabel('$M_r$ [Nm]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%%% magnetotorquer Dipole
figure
plot(T(indexPoint), Dmagn(indexPoint, 1), 'LineWidth', 1); grid on; hold on;
plot(T(indexPoint), Dmagn(indexPoint, 2), 'LineWidth', 1);
plot(T(indexPoint), Dmagn(indexPoint, 3), 'LineWidth', 1);
legend('$D_{x}$', '$D_{y}$', '$D_{z}$')
xlabel('t [s]'); ylabel('D [A$m^2$]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

%% ECLIPSE
indexEcl = find(T >= 2500);
index1 = find(eclipse >=1, 1, 'first');
index2 = find(eclipse >=1, 1, 'last');
%%% is eclipse
figure
plot(T, eclipse, 'LineWidth', 1); grid on; hold on;
xlabel('t [s]'); ylabel('Eclipse [-]')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
ylim([-0.1, 1.1])
set(gca, 'FontSize', 13)

%%% Attitude matrix error
figure
plot(T(indexEcl), squeeze(attErr(1, 1, (indexEcl))), 'LineWidth', 1); grid on; hold on;
plot(T(indexEcl), squeeze(attErr(1, 1, (indexEcl))), 'LineWidth', 1);
plot(T(indexEcl), squeeze(attErr(1, 1, (indexEcl))), 'LineWidth', 1);
plot(T(indexEcl), squeeze(attErr(2, 1, (indexEcl))), 'LineWidth', 1);
plot(T(indexEcl), squeeze(attErr(2, 2, (indexEcl))), 'LineWidth', 1);
plot(T(indexEcl), squeeze(attErr(2, 3, (indexEcl))), 'LineWidth', 1);
plot(T(indexEcl), squeeze(attErr(3, 1, (indexEcl))), 'LineWidth', 1);
plot(T(indexEcl), squeeze(attErr(3, 2, (indexEcl))), 'LineWidth', 1);
plot(T(indexEcl), squeeze(attErr(3, 3, (indexEcl))), 'LineWidth', 1);
xlabel('t [s]'); ylabel('$A_{B/N}$ - $\bar{A}_{B/N}$')
xline(T(index1), '--k', 'start eclipse')
xline(T(index2), '--k', 'end eclipse')
xlow = xlim; xlow = xlow(1); xup = xlim; xup = xup(2);
xticks(floor(linspace(xlow, xup, 6)))
set(gca, 'FontSize', 13)

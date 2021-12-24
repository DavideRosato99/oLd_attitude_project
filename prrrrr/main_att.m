clc;
% controller
C_Kd = 5;
C_Kp = 5;
% attitude dynamics
j0 = diag([0.0109 0.0506 0.0509]);
% reaction wheel
jc = diag([0.00109 0.00109 0.00109]);
R_Km = 0.053;
R_Kd = 0.001;
R_Ka = 0.001;
R_Kn = 0.001;
R_Kl = 0.001;
% Gyroscopes
nd = [0.001 0.001 0.001];
b = [0,0,0]; 
ng = [0.001,0.001,0.001];
D = diag([0.5 1/3 0.25]);
% earth sensor
E_b = [0,0];
nm = [0.001,0.001];
% attitude determination
om_0 = 0.001; % [rad/s]
Kx = 0.012;
Ky = 0.01;
Kz = 0.008;













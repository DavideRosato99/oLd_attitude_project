clear all; close all; clc;
%%
% Ds = diameter of the sun
% Dp = diameter of the planet
% delta_ps = distance between the planet and the sun
% r_eci = spacecraft position vector wrt earth
% d = distance between the center of the umbra and the spacecraft
% s = direction of the sun wrt the earth
% k = penumbra terminator
% csi = umbra terminator
% shadow terminators can be 0,2 or 4:
% case 1: the sc does not enter in the eclpise/shadow zones
% case 2: (penumbra entry(1) and exit(2))
% case 3: (penumbra entry(1), umbra entry(2), umbra exit(3), penumbra exit(4))

%% umbral geometry
X_u = (Dp*delta_ps)/(Ds-dp);
alpha_u = asin(Dp/(2*X_u));

%% penumbral geometry
X_p = (Dp*delta_ps)/(Ds+dp);
alpha_p = asin(Dp/(2*X_u));

%% projection of r_eci on the sun direction
r_s = dot(r_eci,s);

%% distance between the center of the umbra and the spacecraft
d = norm(r_eci) - r_s;
%      ^--- non mi ricordo assolutamente se si faccia cosi la distanza 

%% shadow terminator
% center-penumbra
k = (X_p + abs(r_s))*tan(alpha_p);
% center_umbra 
csi = (X_p - abs(r_s))*tan(alpha_u);

%% bisection method ... jumped
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shadow terminator exists if dot(r_eci,s) < 0
% if abs(d) > k the sc is still in the sunligth
% if csi < abs(d) < k sc is in the penumbra zone
% if abs(d) < csi sc is in the umbra zone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





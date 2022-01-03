function [orb] = car2par(rr, vv, mu)
% car2par - Trasformation from cartesian parameters to Keplerian
%           coordinates
%
% PROTOTYPE
%   [orb]=car2par(rr,vv,mu)
%
% INPUT:
%   rr       double [3x1]   position vector                     [km]
%   vv       double [3x1]   velocity vector                   [km/s]
%   mu       double [1x1]   gravitational parameter       [km^3/s^2]
%
% OUTPUT:
%   orb      double [6x1]   orbital parameters                   [-]
%
% CALLED FUNCTIONS: astroConstants
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

r = norm(rr);
v = norm(vv);
k = [0 0 1]';

a = -mu / (v^2 - 2*mu/r);

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr/r;
e = norm(ee);

i = acos(hh(3)/h);

N = cross(k,hh) / norm(cross(hh,k));

Om = acos(N(1));
if N(2) < 0
    Om = 2*pi - Om;
end

om = acos(dot(N,ee)/e);
if ee(3) < 0
    om = 2*pi - om;
end

if dot(rr,ee)/(r*e) > 1
    theta = 0;
else
    theta = acos(dot(rr,ee)/(r*e));
end

if dot(rr,vv) < 0
    theta = 2*pi - theta;
end

orb = [a e i Om om theta]';


end
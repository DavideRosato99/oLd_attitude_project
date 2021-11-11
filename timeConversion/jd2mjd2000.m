function mjd2000 = jd2mjd2000(jd)

% jd2mjd2000.m - Modified Julian day 2000 number from Julian day number.
%
% PROTOTYPE:
%   mjd2000 = jd2mjd2000(jd)
%
% DESCRIPTION:
%   Returns the modified Julian day 2000 number corresponding to the given
%   Julian day number.
%
% INPUT:
%   jd[1]       Date in Julian Day. The JD (Julian day) count is from 0 at
%               12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%               calendar. The corresponding date in Gregorian calendar is
%               12:00 noon, 24 November -4713.
%
% OUTPUT:
%   mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%               since 01-01-2000, 12:00 noon.
%
% See also mjd20002jd.
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Nicolas Croisard, 16/02/2008, MATLAB, jd2mjd2000.m
%
% CHANGELOG:
%   29/02/2008, REVISION, Camilla Colombo
%   22/04/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% ------------------------- - SpaceART Toolbox - --------------------------


mjd     = jd - 2400000.5;
mjd2000 = mjd - 51544.5;


return
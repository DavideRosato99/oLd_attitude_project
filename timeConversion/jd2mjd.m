function mjd = jd2mjd(jd)

% jd2mjd.m - Modified Julian day number from Julian day number.
%
% PROTOTYPE:
%   mjd = jd2mjd(jd)
%
% DESCRIPTION:
%   Returns the modified Julian day number corresponding to
%   the given Julian day number.
%
% INPUT:
%   jd[1]       Date in Julian Day. The JD (Julian day) count is from 0 at
%               12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%               calendar. The corresponding date in Gregorian calendar is
%               12:00 noon, 24 November -4713.
%
% OUTPUT:
%   mjd[1]      Date in modified Julian Day. The MJD count is from 00:00
%               midnight at the beginning of Wednesday November 17, 1858.
%
% See also mjd2jd.
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Nicolas Croisard, 16/02/2008, MATLAB, jd2mjd.m
%
% CHANGELOG:
%   29/02/2008, REVISION, Camilla Colombo
%   22/04/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------


mjd = jd - 2400000.5;


return
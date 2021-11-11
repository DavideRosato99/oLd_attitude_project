function [hrs, mn, sec] = fracday2hms(fracDay)

% fracday2hms.m - Converts a fraction of day into hours, minutes, and
%	seconds.
%
% PROTPTYPE:
%   [hrs, mn, sec] = fracday2hms(fracDay)
%
% DESCRIPTION:
% 	Converts the fraction of day to hours, minutes, and seconds.
%
% INPUT:
%    fracDay[1] A single real greater or equal to 0 and strictly lower than
%               1.
%
% OUTPUT:
%    hrs[1]     Number of hours as integer greater or equal to 0 and lower
%               or equal to 23.
%    mn[1]      Number of minutes as integer greater or equal to 0 and
%               lower or equal to 59.
%    sec[1]     Number of seconds as a real greater or equal to 0 and
%               strictly lower than 60.
%
% See also hms2fracday.
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%	Nicolas Croisard, 16/02/2008, MATLAB, fracday2hms.m
%
% CHANGELOG:
%   21/02/2008, REVISION, Matteo Ceriotti
%   22/04/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

% Check the input
if nargin ~= 1 || numel(fracDay) ~= 1
    error('FRACDAY2HMS:incorrectInput',...
          'The input should a single element');
end
if fracDay<0 || fracDay>=1
    error('FRACDAY2HMS:incorrectInput',...
          ['The input should be real greater or equal to 0 ,'...
           'and strictly lower than 1']);
end

temp = fracDay*24;
hrs = fix(temp);
mn = fix((temp-hrs)*60);
sec = (temp-hrs-mn/60)*3600;


return
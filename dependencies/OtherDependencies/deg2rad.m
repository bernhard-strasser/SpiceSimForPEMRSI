function rad = deg2rad(deg)
%
% deg2rad Convert from degree to radians.
%
% This function was written by Bernhard Strasser, July 2012.
%
%
%
% Input: 
% -         deg					            ...     An angle in degree.
%
% Output:
% -        rad                              ...     An angle in radians.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None





%% 0. Preparations

% Assign standard values to variables if nothing is passed to function.
if(nargin < 1)
	fprintf('\nError in deg2rad: No Input.\n')
	rad = 0;
	return;
end




%% 1. Convert

rad = pi*deg/180;

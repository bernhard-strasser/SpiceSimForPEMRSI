function [memused,memfree] = memused_linux(quiet_flag)
%
% memused_linux Show memory usage of MATLAB and free memory in GiB.
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function shows the memory that is used by MATLAB in percent of all memory and the free system memory. Should work on all linux systems, probably even on all UNIX
% 
%
%
% memused = memused_linux(quiet_flag)
%
% Input: 
% -         quiet_flag                  ...    If 1, nothing is printed to display.
%
% Output:
% -         memused                     ...    used memory by MATLAB, in percent of all memory
% -         memfree                     ...    free memory on system
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('quiet_flag','var'))
    quiet_flag = 0;
end

% Check if we are on a linux system
[UnixSystem,memfree] = unix('free -g'); UnixSystem = ~UnixSystem;
if(~UnixSystem)
	if(~quiet_flag)
		fprintf('\nWarning: memused_linux cannot run on non-linux systems')
	end
	memused=0;memfree=0;
	return;
end

% 0.3 Definitions
    






%% 1. Find out used memory


% Find out username.
[stat,uname] = unix('whoami');
[stat,uid] = unix('id');
uid = regexp(regexp(uid,'uid=\d+','match'),'\d+','match');
uid = uid{:}; uid = str2num(uid{:});

[stat,memused] = unix(['ps aux | grep "' uname(1:end-1) '\|' num2str(uid) '" | awk ''{print $4"\t"$11}'' | grep -i "matlab" | cut -f 1']);     % unix just performs the unix-command.
memused = sum(str2num(memused));                                                               % str2double does not work

% free changed from version 3.3.9 to 3.3.10 outputting different fields
[stat,freeVersion] = unix('free -V');
freeVersion = freeVersion(regexp(freeVersion,'[\d+\.]+\d+'):end);
freeVersion = strsplit(freeVersion,'.');
freeVersion = str2double(regexprep(freeVersion,'\D',''));

if(any(isnan(freeVersion)) || freeVersion(1) > 3 || freeVersion(2) > 3 || freeVersion(3) > 9)
    [stat,memfree] = unix('free -m | awk ''NR==2'' | awk ''{print $7}''');     
else
    [stat,memfree] = unix('free -m | awk ''NR==3'' | awk ''{print $4}''');  % unix just performs the unix-command. perform free -g command,
end                                                                         % take 3rd line of that, and 4th column, i.e. "-/+ buffers/cache" of free column
clear stat
memfree = sum(str2num(memfree));                                                               % str2double does not work


%% 2. Display used memory

if(~quiet_flag)
    display([ char(10) num2str(memused) '% of total memory is used by MATLAB.' char(10)])
    display([ char(10) num2str(memfree) ' MB of memory available.' char(10)])    
end




%% 3. Postparations

% fclose(fid)







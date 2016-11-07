function vtbchk
%VTBCHK Checks for updates on a regular basis.
% chkskip tells how often to run check.

%Joseph C. Slater, April 2008
chkskip=7;% number of days to go without checking again.
curpath=pwd;
vtbdir=which('vtb1_1.m');vtbdir=vtbdir(1:(length(vtbdir)-8));
cd(vtbdir)
if exist('chkdatestamp.txt')==0
	chkdatestamp=0;
else
	[chkdatestamp,status]=urlread(['file:///' fullfile(vtbdir,'chkdatestamp.txt')]);
end


if (str2double(chkdatestamp)<(now-chkskip))
	[insdatestamp,status]=urlread(['file:///' fullfile(vtbdir,'vtbdatestamp.txt')]);
	[curdatestamp,status]=urlread('http://www.cs.wright.edu/people/faculty/jslater/vtoolbox/vtoolbox/vtbdatestamp.txt');
	if (str2double(insdatestamp)<str2double(curdatestamp))
	vtbud
	disp('Run vtbud at any time you are online to check for updates to the Engineering Vibration Toolbox.')
    end
    if status==0
    disp('Engineering Vibration Toolbox update checking not working. No network connection.')
    disp('Run ''vtbchk'' while on online to check for updates.')
    disp(['Automatic check again in ' str2num(chkskip) ' days.'])
    end
end

	

fid = fopen('chkdatestamp.txt','wt');
fprintf(fid,'%s',num2str(now));
fclose(fid);
cd(curpath)


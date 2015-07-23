function vtbud
% Displays updates to Engineering Vibration Toolbox after 1/1/98
% Launches browser to show current updates available online
%
% Prompts for automated upating and installation
%
% If Matlab is not set up to launch your browser, connect to
% web https://raw.githubusercontent.com/vibrationtoolbox/vtoolbox/master/vtbud.txt
% manually.
%

% type vtbud.txt
%          Copyright Joseph C. Slater, October 2006

head = 'http://vibrationtoolbox.com/';
ziploc = 'https://github.com/vibrationtoolbox/vtoolbox/archive/master.zip';
webpageloc= 'http://vibrationtoolbox.github.io';

if ~strcmp(java.lang.System.getProperty( 'java.awt.headless' ),'true')


vtbdir=which('vtb1_1.m');
curdir=pwd;
license='No';
if isempty(vtbdir)
    answer=questdlg('Do you want to install the Engineering Vibration Toolbox?');
    if strcmp(answer,'Yes')
        license=['LICENSE The Engineering Vibration Toolbox is licensed free of charge for educational use '... 
                'For professional use, users should contact the Engineering Vibration Toolbox' ... 
                'author directly. Joseph C. Slater is the copyright holder of the Engineering '... 
                'Vibration Toolbox. Neither the author, Prentice Hall, nor Wright State University '... 
                'make any warranty with regard to merchantability or fitness for any given '...
                'purpose with regard to the software. All rights are retained. No '...
                'permission is given to anyone other than myself, the MathWorks and '...
                'Prentice Hall to distribute this software in any manner whatsoever. '...
                'Usage without ownership of Engineering Vibration by Daniel Inman is prohibited without '...
                'the express permission of the author of the software.'];
        license=questdlg(license,'Do you agree?','No');
    end

    if strcmp(answer,'Yes')&&strcmp(license,'Yes')
        vtbupdir = uigetdir(pwd,'Choose a directory to install the Engineering Vibration Toolbox in.') ;
        if ~isempty(findstr(lower(vtbupdir),'matlab'))
            ans1=questdlg('Installation of the Engineering Vibration Toolbox inside the Matlab application directory can cuase it to be deleted during future matlab updates. Are you sure you want to install here?','Are you sure?','No');
            if strcmp(ans1,'No')
                vtbupdir = uigetdir(pwd,'Choose a directory to install the Engineering Vibration Toolbox in.') ;
            end
        end
    end
else
    vtbdir=vtbdir(1:(length(vtbdir)-8));
    vtbupdir=vtbdir(1:length(vtbdir)-9);
%    web vtbud.txt
    
    insdatestamp=urlread(['file:///' fullfile(vtbdir,'vtbdatestamp.txt')]);
    curdatestamp=urlread([head 'chkdatestamp.txt']);
    if str2num(insdatestamp)<str2num(curdatestamp)&~strcmp(java.lang.System.getProperty( 'java.awt.headless' ),'true')
    web([head 'vtbud.txt']);
    answer=questdlg('You do not have the most recent version of the Engineering Vibration Toolbox. Please, review the recent updates to determine if you want to update. Do you want to update the Engineering Vibration Toolbox?','Update now?','Update','No','Cancel','tex');
        if strcmp(answer,'Update')
            answer='Yes';
            update='Yes';
            license=['LICENSE The Engineering Vibration Toolbox is licensed free of charge for educational use '... 
                    'For professional use, users should contact the Engineering Vibration Toolbox' ... 
                    'author directly. Joseph C. Slater is the copyright holder of the Engineering '... 
                    'Vibration Toolbox. Neither the author, Prentice Hall, nor Wright State University '... 
                    'make any warranty with regard to merchantability or fitness for any given '...
                    'purpose with regard to the software. All rights are retained. No '...
                    'permission is given to anyone other than myself, the MathWorks and '...
                    'Prentice Hall to distribute this software in any manner whatsoever. '...
                    'Usage without ownership of Engineering Vibration by Daniel Inman is prohibited without '...
                    'the express permission of the author of the software.'];
            license=questdlg(license,'Do you agree?','No');

        elseif strcmp(answer,'No')
            answer=questdlg('Do you want to remove the Engineering Vibration Toolbox?','Remove?','No');
            if strcmp(answer,'Yes')
                answer='Cancel';
                rmdir('vtoolbox','s')
                rmpath([vtbupdir 'vtoolbox'])     
            end
        end
    else
        answer=questdlg('You have the most recent version of the Engineering Vibration Toolbox.  Do you want to remove the Engineering Vibration Toolbox?','Remove?','No');
        
        if strcmp(answer,'Yes')
            annswer=questdlg('Are you sure you want to remove the Engineering Vibration Toolbox?','Really Remove?','No');
            if strcmp(annswer,'Yes')
                answer='Cancel';
                rmdir('vtoolbox','s')
                rmpath([vtbupdir 'vtoolbox'])    
            end
        end        
    end    
end
if strcmp(answer,'No')
    answer='Cancel';
end

if strcmp(answer,'Yes')&&strcmp(license,'Yes')
    cd(vtbupdir)
    if exist('vtoolbox.zip','file')==2
        answer=questdlg('Do you want to install from the internet or vtoolbox.zip on your hard drive?',...
        'Install from?','Hard drive','Internet','Internet');
        if strcmp(answer,'Internet')
            delete('vtoolbox.zip')
            urlwrite(ziploc,'vtoolbox.zip')
            unzip('vtoolbox.zip')
            delete('vtoolbox.zip')
        else
            unzip('vtoolbox.zip',pwd)
        end
    else
        [weby,st]=urlread(webpageloc);
        if st==0
            msgdlg('Without a downloaded vtoolbox.zip file or an internet connection the installation cannot continue.','Error')
            return
        end
        urlwrite(ziploc,'vtoolbox.zip')
        unzip('vtoolbox.zip')
        delete('vtoolbox.zip')
    end
cd('vtoolbox');
addpath(pwd,'-end');
vtbdir=pwd;
aa=savepath;
if aa==1
    if isempty(findstr('vtoolbox',path))
    if isunix&&~strcmp(computer,'MACI')&&~strcmp(update,'Yes')
        astartmod=questdlg(['Path not saved. You will need to add the line ''addpath(''' pwd ''')'' to your startup.m file. Do you want me to attempt to do this?']) ;
        if strcmp(astartmod,'Yes')
            ucommand=['!echo addpath\(\''' pwd '\''' ',' '\''' '-' 'end' '\''' '\) >> ~/startup.m'];
            eval(ucommand);
            msgbox('The last line of the file startup.m in your home directory should now be set to add vtoolbox to your path each time you run Matlab. In order for the Engineering Vibration Toolbox to work, you must always run Matlab from your home directory. ');
        else
            msgbox('You must add the command ''''addpath ''vtoolbox'' -end'''' to your startup.m file in order to use the Engineering Vibration Toolbox.');
        end
        
    else
        warndlg('Path not saved. You will have to manually set your startup path to include the vtoolbox directory.');
        pathtool
    end
    end
end
cd(curdir)
msgbox(['Installation of the Engineering Vibration Toolbox in ' vtbdir ' is complete. Please delete all copies of vtbud.m that are not in the vtoolbox directory.'],'Install Complete');
end
if strcmp(license,'No')&&~(strcmp(answer,'Cancel')||strcmp(answer,'Remove'))
    msgbox('Please note that usage without agreement to the license is prohibited.','Usage prohibited.');
else
    if ~strcmp(answer,'Cancel')
    doc vtoolbox
    end
end
else
    disp('vtbud only works when you have a GUI running. Sorry.')
end


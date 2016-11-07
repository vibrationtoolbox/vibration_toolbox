%Simple script to add the vtoolbox directory to the matlab path


vtbupdir = uigetdir(pwd,'Choose the Engineering Vibration Toolbox (vtoolbox) directory.'); 
cd(vtbupdir)
addpath(pwd,'-end')
aa=savepath;
if aa==1
    if isunix&&~strcmp(computer,'MACI')
        astartmod=questdlg(['Path not saved. You will need to add the line ''addpath(''' pwd ''')'' to your startup.m file. Do you want me to attempt to do this?']) ;
        if strcmp(astartmod,'Yes')
            ucommand=['!echo addpath\(\''' pwd '\''' ',' '\''' '-' 'end' '\''' '\) >> ~/startup.m'];
            eval(ucommand);
            msgbox('The last line of the file startup.m in your home directory should now be set to add vtoolbox to your path each time you run Matlab. In order for the Engineering Vibration Toolbox to work, you must always run Matlab from your home directory. ')
        end
    else
        warndlg('Path not saved. You will have to manually set your startup path to include the vtoolbox directory.')
        pathtool
    end

end

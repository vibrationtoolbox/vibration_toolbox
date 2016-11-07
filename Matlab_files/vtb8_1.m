function vtb8_1
% VTB8_1  [node,ncon,zero,force]=VTB8_1
% Makes input file for VTB8_2
% I suggest you print the file vtb8read.txt

%	Joseph C. Slater, 6-17-90
%	Copyright (c) 1990-94 by Joseph C. Slater

%       11 Nov 94 -- fixed graphics to work with Matlab 4.0

clc
home
aa=version;ll=length(aa);
figure
grid on
loc=input('Enter x and y location of node. (eg.: [x y]) ');
node(1,:)=loc;
plot(node(:,1),node(:,2),'*b')
axis('square')
nnum=['1  '];
for i=2:1000
  loc=input('Enter x and y location of node (ie. [x y]) or 0 to end. ');
  if loc==0 & length(loc)==1 ,break,end
  node(i,:)=loc;
  len=max([max(node(:,1))-min(node(:,1)) max(node(:,2))-min(node(:,2))]);
%  xl=(max(node(:,1))+min(node(:,1)))/2-.7*length;
%  xh=(max(node(:,1))+min(node(:,1)))/2+.7*length;
%  yl=(max(node(:,2))+min(node(:,2)))/2-.7*length;
%  yh=(max(node(:,2))+min(node(:,2)))/2+.7*length;
%%%
  d=node;
  xlo=min([node(:,1);d(:,1)]);
  xho=max([node(:,1);d(:,1)]);
  ylo=min([node(:,2);d(:,2)]);
  yho=max([node(:,2);d(:,2)]);
  xsp=xho-xlo;
  ysp=yho-ylo;
  if .618*xsp > ysp
      xh=xho+.1*xsp;
      xl=xlo-.1*xsp;
	  yh=(yho+ylo)/2+(xh-xl)*.618/2;
	  yl=(yho+ylo)/2-(xh-xl)*.618/2;
	else
      yh=yho+.1*ysp;
      yl=ylo-.1*ysp;
      xh=(xho+xlo)/2+.5*(yh-yl)/.618;
      xl=(xho+xlo)/2-.5*(yh-yl)/.618;
  end
%%%
  dx=.01*(xh-xl);
  dy=.035*(yh-yl);
  plot(node(:,1),node(:,2),'*b')
  j=1:i;
  if i < 10
      istr=[num2str(i) '  '];
	elseif i < 100
	  istr=[num2str(i) ' '];
	else
	  istr=[num2str(i)];
  end	
  nnum=[nnum;istr];
  text(node(j,1)'+dx,node(j,2)'+dy,nnum)
  axis([xl xh yl yh])
  axis('image')
end
hold on
clc
%disp('Do you have a pointing device such as a mouse or trackball? (y/n)')
%disp('(Arrow keys may be sufficient)')
%point=input(' ','s');
point='y';
%connecting nodes section
clc
home
  disp(' Pick nodes to connect with elements.')
if point=='y'
  disp(' (Use pointing device or arrow keys and return)')
 else
  disp(' (Enter node numbers one at a time)')
end
  pause(1)
answer2='n';
for i=1:1000
%  clc
%  home
if i~=1
  disp(' Pick nodes to connect with elements.')
end
  if point=='y'
    [x1 y1]=ginput(1);
    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
    [dsq,nodenum1]=min(dis);
   else
    clc,home
    nodenum1=input(' Enter node number 1: ');
  end
  plot(node(nodenum1,1),node(nodenum1,2),'*w')
  plot(node(nodenum1,1),node(nodenum1,2),'or')
  if point=='y'
    [x2 y2]=ginput(1);
    dis=(node(:,1)-x2).^2+(node(:,2)-y2).^2;
    [dsq,nodenum2]=min(dis);
   else
    clc,home
    nodenum2=input(' Enter node number 2: ');
  end
  plot(node(nodenum2,1),node(nodenum2,2),'*w')
  plot(node(nodenum2,1),node(nodenum2,2),'or')
  pause(.04)
  plot([node(nodenum1,1) node(nodenum2,1)],[node(nodenum1,2) node(nodenum2,2)],'-b')
  plot(node(nodenum1,1),node(nodenum1,2),'ow')
  plot(node(nodenum1,1),node(nodenum1,2),'*b')
  plot(node(nodenum2,1),node(nodenum2,2),'ow')
  plot(node(nodenum2,1),node(nodenum2,2),'*b')
  clc
  home

  if i>1
  answer2=input(' Same properties as previous element? (y/n) ','s');
  end
  if answer2=='n'
    clc
    E=input(' Enter the modulus of elasticity of the member. ');
    G=input(' Enter the shear modulus of the member (zero for EB beam). ');
    I=input(' Enter the moment of area of the member. ');
    A=input(' Enter the cross sectional area of the member. ');
    Rho=input(' Enter the density per unit length of the member. ');
  end
  ncon(i,:)=[nodenum1 nodenum2 E A I G Rho];
  answer=input(' Enter another element? (y/n) ','s');
  if answer~='y',break,end
end



% adding concentrated masses and inertias
conm=[];
for i=1:1000
clc
home
  if i==1
    answer=input(' Add concentrated masses and rotational inertias? (y/n) ','s');
  end
  if i>1
    answer=input(' Add more concentrated masses and rotational inertias? (y/n) ','s');
  end
  if answer~='y',break,end
  disp(' ')
  disp(' Pick node to add mass/rotational inertia to.')
  pause(.5)
  if point=='y'
    [x1 y1]=ginput(1);
    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
    [dsq,nodenum]=min(dis);
   else
    nodenum=input(' Enter node number: ');
  end
  plot(node(nodenum,1),node(nodenum,2),'*w')
  plot(node(nodenum,1),node(nodenum,2),'xr')
  answer=input(' Add mass or rotational inertia?(m,i,n(one)) ','s');
  massval=input(' Enter magnitude of mass/inertia. ');
  conm(i,:)=[0 0 0];
  conm(i,1)=nodenum;

  if answer=='m'
    conm(i,2)=massval;
  end
  if answer=='i'
    conm(i,3)=massval;
  end

  plot(node(nodenum,1),node(nodenum,2),'xi')
  plot(node(nodenum,1),node(nodenum,2),'*b')
end


% zeroing of displacements
zero=[];ij=0;
for i=1:1000
clc
home
  if i==1
    answer=input(' Add boundary conditions? (y/n) ','s');
  end
  if i>1
    answer=input(' Zero another displacement? (y/n) ','s');
  end
  if answer~='y',break,end
  disp(' ')
  disp(' Pick node to zero')
  pause(1.)
  if point=='y'
    [x1 y1]=ginput(1);
    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
    [dsq,nodenum]=min(dis);
   else
    nodenum=input(' Enter node number: ');
  end
  plot(node(nodenum,1),node(nodenum,2),'*w')
  plot(node(nodenum,1),node(nodenum,2),'xr')
  disp(' ')
  disp(' ')
  disp(' Zero which displacement(s)?(x,y,t(heta),n(one))')
  disp(' (ie xt for x and theta)');
  answern=input(' ','s');
  sanswern=size(answern);
for ii=1:sanswern(2)
    ij=ij+1;
    answer=answern(ii);
    if answer=='x'
%       plot(node(nodenum,1),node(nodenum,2),'xi')
       plot(node(nodenum,1),node(nodenum,2),'*b')
       zero(ij,:)=[nodenum 1];
    end
    if answer=='y'
%       plot(node(nodenum,1),node(nodenum,2),'xi')
       plot(node(nodenum,1),node(nodenum,2),'*b')
       zero(ij,:)=[nodenum 2];
    end
    if answer=='t'
%       plot(node(nodenum,1),node(nodenum,2),'xi')
       plot(node(nodenum,1),node(nodenum,2),'*b')
       zero(ij,:)=[nodenum 3];
    end
    if answer=='n'
%       plot(node(nodenum,1),node(nodenum,2),'xi')
       plot(node(nodenum,1),node(nodenum,2),'*b')
    end
  end
end
force=[];
% adding loads
for i=1:1000
clc
home
  answer=input(' Add loads? (y/n) ','s');
  if answer~='y',break,end
  disp(' ')
  disp(' Pick node to load')
  pause(1.5)
  if point=='y'
    [x1 y1]=ginput(1);
    dis=(node(:,1)-x1).^2+(node(:,2)-y1).^2;
    [dsq,nodenum]=min(dis);
   else
    nodenum=input(' Enter node number: ');
  end
  plot(node(nodenum,1),node(nodenum,2),'*w')
  plot(node(nodenum,1),node(nodenum,2),'xr')
  answer=input('Load which displacement?(x,y,t(heta),n(one)) ','s');

  if answer=='x'
    loadval=input(' Enter magnitude of load. ');
    force(i,:)=[nodenum 1 loadval];
  end

  if answer=='y'
    loadval=input(' Enter magnitude of load. ');
    force(i,:)=[nodenum 2 loadval];
  end

  if answer=='theta'
    loadval=input(' Enter magnitude of load. ');
    force(i,:)=[nodenum 3 loadval];
  end

  answer=input(' Load another node? (y/n) ','s');
  if answer~='y',break,end
%  plot(node(nodenum,1),node(nodenum,2),'xi')
  plot(node(nodenum,1),node(nodenum,2),'*b')
end
node;
ncon;
zero;
force;

%path(path,pwd)
answer=input('Save configuration file (Else all will have been in vain)? (y/n) ','s');
%if answer=='y'
%  filename=input(' Enter name of configuration file. ','s');
%  eval(['save ',filename,'.con',' node',' ncon',' zero',' force',' conm']);
%  answer=input(' Run analysis? (y/n) ','s');
%  if answer=='y'
%    vtb8_2(filename);
%    vtb8_2(node,ncon,zero,force,conm);
%  end
%end
if answer=='y'
  [filename,pathname]=uiputfile('projectname.con','Save as:');
  sfilename=size(filename);
  path(path,pathname)
  lfilename=sfilename(2);
  if filename(lfilename-3)=='.'
    filename=[filename(1:lfilename-4)  '.con'];
    projectname=filename(1:lfilename-4);
   else
    filename=[filename '.con'];
    projectname=filename(1:lfilename-4);
  end
  pathname;
  sizepath=size(pathname);
  shortpathname=pathname(1:sizepath(2)-1);
  lsp=size(findstr(shortpathname,':'));
  if strcmp(computer,'MAC2') & lsp(1)==0
    shortpathname=[shortpathname ':'];
  end
  %path(shortpathname)
  cdpath=['cd ' '''' shortpathname '''' ];% Crazy quotes allow spaces
  %                                       % in directory names.
  eval(cdpath)
  eval(['save ',filename,' node',' ncon',' zero',' force',' conm']);
  answer=input(' Run analysis? (y/n) ','s');
  if answer=='y'
    vtb8_2(projectname);
  end
end

%Automatically check for updates
vtbchk

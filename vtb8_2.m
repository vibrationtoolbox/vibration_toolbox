function [xout,fout]=VTB8_2(node,ncon,zero,force,conm)
%VTB8_2 [x,f]=VTB8_2(node,ncon,zero,force,conm)
%       [x,f]=vtb8_2('filename');
%    Timoshenko 2-D beam finite element code.
%
%    node=[x1 y1;x2 y2;...]
%    ncon=[node1 node2 E A I G Rho;...]
%        Where 'node1' and 'node2' are connected by an element,
%        'E' is Young's modulus, 'A' is the cross sectional area,
%        'I' is the moment of area, 'G' is the shear modulus
%        and 'Rho' is the density per unit length
%        (set 'G' equal to zero to ignore shear deformation).
%        For pure truss elements, set I=0 and zero all rotations.
%    zero=[node# dof#;...]
%        'dof#' is the degree of freedom at node 'node#'
%            to constrain or load.
%        'dof#' numbers [1 2 3] correspond to [x y theta]
%            force is the magnitude of the load.
%    force=[node# dof# force]
%    conm=[node# mass rotational inertia]
%        Where 'node#' is the node at which the mass of magnitude
%        'mass' is located.
%    All rotations are positive counter clockwise.
%    Four methods exist for creating this data.
%      1) Use the program vtb8_1.
%      2) Type clear.  Enter the data interactively.
%         Save a file with the file extension .con.
%         i.e. type: save beam1.con
%         Type: beam
%      3) Type clear.  Enter the data interactively.
%         Type: [x,f]=vtb8_2(node,ncon,zero,force)
%      4) Create a script 'm' file including the definitions.
%         Add the line:
%         save 'filename.con'
%         to the end of the file.
%         Execute the script file.  
%         Type: [x,f]=vtb8_2(node,ncon,zero,force,conm);
%         or
%         [x,f]=vtb8_2('filename'); % Note, no extension.
%    An example file named vtb8_e1.m is available.
%
%    I suggest you print the file vtb8read.txt.
clc
home
%help readme10
%disp(' ')
%disp(' Hit return to continue')
%pause

clc,home

if nargin==0
  [filename,pathname]=uigetfile('*.con','Open Connectivity File');
  sizepath=size(pathname);
  shortpathname=pathname(1:sizepath(2)-1);
  lsp=size(findstr(shortpathname,':'));
  if strcmp(computer,'MAC2') & lsp(1)==0
    shortpathname=[shortpathname ':'];
  end
  cdpath=['cd ' '''' shortpathname '''' ];% Crazy quotes allow spaces
  %                                       % in directory names.
  eval(cdpath)
  sfilename=size(filename);
  path(path,pathname);
  oflag=1;
  lfilename=sfilename(2);
  if filename(lfilename-3)=='.'
    filename=[filename(1:lfilename-4)  '.con'];
	projectname=filename(1:lfilename-4);
   else
    filename=[filename '.con'];
	projectname=filename(1:lfilename-4);
  end
  eval(['load ',filename, ' -mat']);
end

if nargin==1
  disp(' Loading configuration data from file.')
  projectname=node;
  eval(['load ',projectname,'.con -mat']);
end

clc
home
disp(' Constructing Global Mass and Stiffness Matrices')
disp(' ')
snode=size(node);
k=zeros(3*snode(1));
m=zeros(3*snode(1));
sncon=size(ncon);
szero=size(zero);
% Assembly of mass and stiffness matrices.

for ii=1:sncon(1)
  iis=num2str(ii);
  ke=zeros(6,6);
  me=zeros(6,6);
  n1=ncon(ii,1);
  n2=ncon(ii,2);
  x1=node(n1,1);
  y1=node(n1,2);
  x2=node(n2,1);
  y2=node(n2,2);
  theta=atan2(y2-y1,x2-x1);
  s=sin(theta);
  c=cos(theta);
  le=sqrt((y2-y1)^2+(x2-x1)^2);
  E=ncon(ii,3);
  A=ncon(ii,4);
  I=ncon(ii,5);
  G=ncon(ii,6);
  Rho=ncon(ii,7);
  if I==0
    I=1e-8*A;
  end
  R=A*le^2/I;
  alpha=1.5;

  %  alpha is equal to the ratio of the maximum to the average
  %    shear stress for the assumed stress distribution 
  %    through the depth of the beam.  3/2 is the value for 
  %    a rectangular cross section.  
theta;
  S=G*A*le^2/12/E/I/alpha;
  if S<1e-10
    mes1=[' No Shear Stiffness in Element ',iis,'.'];
    disp(mes1);
    ke(1,1)=R*c^2+12*s^2;
    ke(1,2)=c*s*(R-12);
    ke(2,2)=R*s^2+12*c^2;
    ke(1,3)=-6*le*s;
    ke(2,3)=6*le*c;
    ke(3,3)=4*le^2;
    ke(3,6)=2*le^2;
  else
    ke(1,1)=R*c^2+12/(1+1/S)*s^2;
    ke(1,2)=c*s*(R-12/(1+1/S));
    ke(2,2)=R*s^2+12/(1+1/S)*c^2;
    ke(1,3)=-6*le*s/(1+1/S);
    ke(2,3)=6*le*c/(1+1/S);
    ke(3,3)=4*le^2*(1+1/4/S)/(1+1/S);
    ke(3,6)=2*le^2*(1-1/2/S)/(1+1/S);
  end
  ke(1,4)=-ke(1,1);
  ke(1,5)=-ke(1,2);
  ke(1,6)=ke(1,3);
  ke(2,4)=-ke(1,2);
  ke(2,5)=-ke(2,2);
  ke(2,6)=ke(2,3);
  ke(3,4)=-ke(1,3);
  ke(3,5)=-ke(2,3);
  ke(4,4:5)=ke(1,1:2);
  ke(4,6)=-ke(1,3);
  ke(5,5)=ke(2,2);
  ke(5,6)=-ke(2,6);
  ke(6,6)=ke(3,3);
  ke=ke*E*I/le^3;
  ke=ke+ke'-diag(diag(ke));

  me(1,1)=140*c^2+156*s^2;
  me(1,2)=-16*c*s;
  me(1,3)=-22*le*s;
  me(1,4)=70*c^2+54*s^2;
  me(1,5)=16*c*s;
  me(1,6)=13*le*s;
  me(2,2)=140*s^2+156*c^2;
  me(2,3)=22*le*c;
  me(2,4)=me(1,5);
  me(2,5)=70*s^2+54*c^2;
  me(2,6)=-13*le*c;
  me(3,3)=4*le^2;
  me(3,4)=-me(1,6);
  me(3,5)=-me(2,6);
  me(3,6)=-3*le^2;
  me(4,4)=me(1,1);
  me(4,5)=me(1,2);
  me(4,6)=-me(1,3);
  me(5,5)=me(2,2);
  me(5,6)=-me(2,3);
  me(6,6)=me(3,3);

  me=me*Rho*le/420;
  me=me+me'-diag(diag(me));
theta;
ke;
me;

  p1=(n1-1)*3+1;
  p2=(n2-1)*3+1;
  k(p1:p1+2,p1:p1+2)=k(p1:p1+2,p1:p1+2)+ke(1:3,1:3);
  k(p1:p1+2,p2:p2+2)=k(p1:p1+2,p2:p2+2)+ke(1:3,4:6);
  k(p2:p2+2,p1:p1+2)=k(p2:p2+2,p1:p1+2)+ke(4:6,1:3);
  k(p2:p2+2,p2:p2+2)=k(p2:p2+2,p2:p2+2)+ke(4:6,4:6);

  m(p1:p1+2,p1:p1+2)=m(p1:p1+2,p1:p1+2)+me(1:3,1:3);
  m(p1:p1+2,p2:p2+2)=m(p1:p1+2,p2:p2+2)+me(1:3,4:6);
  m(p2:p2+2,p1:p1+2)=m(p2:p2+2,p1:p1+2)+me(4:6,1:3);
  m(p2:p2+2,p2:p2+2)=m(p2:p2+2,p2:p2+2)+me(4:6,4:6);


end

% Adding concentrated masses.

sconm=size(conm);
if sconm(1)~=0
  sc=size(conm);
  if sc(2) ==2
    conm(1,3)=0;
  end
  for i=1:sconm(1)
    loc1=(conm(i,1)-1)*3+1;
    loc2=(conm(i,1)-1)*3+2;
    loc3=(conm(i,1)-1)*3+3;
    m(loc1,loc1)=conm(i,2)+m(loc1,loc1);
    m(loc2,loc2)=conm(i,2)+m(loc2,loc2);
    m(loc3,loc3)=conm(i,3)+m(loc3,loc3);
  end
end

% Zeroing stiffness matrix and mass matrix.
clc
home
disp(' Applying Boundary Conditions.')
disp(' ')
k1=k;
m1=m;
if length(zero)~=0
  np=(zero(:,1)-1)*3+zero(:,2);
  np=sort(np);
  p=1:3*snode;
  p(np)=p(np)*0;
  p=sort(p);
  p=p(length(np)+1:length(p));
  k1=k(p,p);
  m1=m(p,p);
end
clc,home
selx=1:3:snode(1)*3;
sely=2:3:snode(1)*3;
selt=3:3:snode(1)*3;
answer=input('Do you want to do a static or dynamic analysis? (s/d) ','s');
if answer=='s' | answer=='S'
  pf=(force(:,1)-1)*3+force(:,2);
  f=zeros(snode(1)*3,1);
  f(pf)=force(:,3);
  f1=zeros(length(k1),1);
  f1=f(p);
  x1=k1\f1;
  x=zeros(snode(1)*3,1);
  x(p)=x1;
  f=k*x;
  xx=[x(selx) x(sely) x(selt)];
  ff=[f(selx) f(sely) f(selt)];
else
  [vk,dk]=eig(k1);
  flag=0;
  for ij=1:length(k1)
    clc,disp('Working')
    if dk(ij,ij)<1e-14
      flag=1;
      keig=num2str(dk(ij,ij));
      dk(ij,ij)=1e-14;
      disp(' Numerical roundoff error occurred')
      disp(' Eigenvalue of stiffness matrix changed from')
      disp([' ' keig ' to 0'])
    end
  end

  if flag == 1
     disp(' '),disp(' Press return to continue'),pause
     disp(' Rebuilding corrected stiffness matrix'),pause(1)
     k1=vk*dk*vk';
  end
  [x1,wsq]=eig(k1,m1);
  f1=wsq.^.5;
  [fs1,sf1]=sort(diag(f1));
  f1=diag(fs1);
  x1=x1(:,sf1);
  snode;
  x=zeros(snode(1)*3,snode(1)*3-szero(1));
  if length(zero)~=0
    for ic=1:length(x1)
      x(p,ic)=x1(:,ic);
    end
  else
    for ic=1:length(x1)
      x(:,ic)=x1(:,ic);
    end
  end
  ff=diag(f1);
  f=ff;
  xx=x;
end

clc,home

answer1=input('Save Results(y)? (y/n) ','s');
if answer1~='n'
  if exist('pathname')
    sizepath=size(pathname);
    shortpathname=pathname(1:sizepath(2)-1);
    lsp=size(findstr(shortpathname,':'));
    if strcmp(computer,'MAC2') & lsp(1)==0
      shortpathname=[shortpathname ':'];
    end
    cdpath=['cd ' '''' shortpathname '''' ];% Crazy quotes allow spaces
    %                                       % in directory names.
    eval(cdpath)
  end
  if exist('projectname')==0
    [filename,pathname]=uiputfile('projectname.out','Save as:');
    sfilename=size(filename);
    path(path,pathname)
    lfilename=sfilename(2);
    if filename(lfilename-3)=='.'
      filename=[filename(1:lfilename-4)  '.out'];
      projectname=filename(1:lfilename-4);
     else
      filename=[filename '.out'];
	  projectname=filename(1:lfilename-4);
    end
  end
  eval(['save ',projectname,'.out',' x',' f']);
  disp('Saving data')
end

answer2=input('Save Equations(y)? (y/n) ','s');
if answer2~='n'
  if exist('pathname')
    sizepath=size(pathname);
    shortpathname=pathname(1:sizepath(2)-1);
    lsp=size(findstr(shortpathname,':'));
    if strcmp(computer,'MAC2') & lsp(1)==0
      shortpathname=[shortpathname ':'];
    end
    cdpath=['cd ' '''' shortpathname '''' ];% Crazy quotes allow spaces
    %                                       % in directory names.
    eval(cdpath)
  end
  if exist('projectname')==0
    [filename,pathname]=uiputfile('projectname.eqn','Save as:');
    sfilename=size(filename);
    path(path,pathname)
    lfilename=sfilename(2);
    if filename(lfilename-3)=='.'
      filename=[filename(1:lfilename-4)  '.eqn'];
      projectname=filename(1:lfilename-4);
     else
      filename=[filename '.eqn'];
	  projectname=filename(1:lfilename-4);
    end
  end
  p=p';
  eval(['save ',projectname,'.eqn',' k1',' m1',' x1',' f1',' p']);
  disp('Saving data')
end

answer=input('Show results graphically(y)? (y/n) ','s');
if answer~='n'
  vtb8_3(node,x,zero,ncon,p,f);
end
nargout;
if nargout==0
 return % Suppress Output
end

xout=x;
fout=f;

%Automatically check for updates
vtbchk

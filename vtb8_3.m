function vtb8_3(node,x,zero,ncon,p,f)
%VTB8_3 allows you to view the results of VTB8_2
% I suggest you print the file vtb8read.

if nargin==0
  [filename,pathname]=uigetfile('*.out','Open Project Output File');
  sfilename=size(filename);
  sizepath=size(pathname);
  shortpathname=pathname(1:sizepath(2)-1);
  lsp=size(findstr(shortpathname,':'));
  if strcmp(computer,'MAC2') & lsp(1)==0
    shortpathname=[shortpathname ':'];
  end
  cdpath=['cd ' '''' shortpathname '''' ];% Crazy quotes allow spaces
  %                                       % in directory names.
  eval(cdpath)
  oflag=1;
  lfilename=sfilename(2);
  if filename(lfilename-3)=='.'
    filename=[filename(1:lfilename-4)  '.con'];
	projectname=filename(1:lfilename-4);
   else
    filename=[filename '.con'];
	projectname=filename(1:lfilename-4);
  end
  eval(['load ',projectname,'.con -mat']);
  eval(['load ',projectname,'.eqn -mat']);
  eval(['load ',projectname,'.out -mat']);
end
answer='c';
snode=size(node);
nnum=[];
flag3=0;
for i=1:snode(1)
 stri=num2str(i);
 if length(stri)==1
  stri=[stri ' '];
 end
 if length(stri)==2
  stri=[stri ' '];
 end
 nnum=[nnum;stri];
end
figure
ndof=size(node);
ndof=ndof(1)*3;
mn=1;
xx=x;
%hold off
scale=1;
nmodes=size(xx);
nmodes=nmodes(2);
if nmodes>1
   flag=1;
  else
   flag=0;
end
for ii=1:1000
  for jj=1:1000
    clc
    disp('')
    if flag==1
      msg1=['Present mode number is: ' num2str(mn) ' of ' num2str(nmodes)];
      disp(msg1)
    end
    msg2=['Present natural frequency is: ' num2str(f(mn)/2/pi) ' Hz (' num2str(f(mn)) ' rad/s)'];
    disp(msg2)
    msg2=['Present scale factor is: ' num2str(scale)];
    disp(msg2)
    disp('')
    disp('')

    disp('Select settings.')
    if flag==0
      disp('Only one deformation found.')
      disp('Selection a is disabled.')
    end
	
    disp('')
    disp('     a) Select Mode Number')
    disp('     b) Scale Deformation')
    disp('     c) Show Deformation')
	disp('     d) Print Deformation')
    disp('     e) Exit')
    answer=input('Enter Choice: ','s');
    if answer=='e' | answer=='c'
      break
    end
    if answer=='a' & flag==1
      clc
      mnold=mn;
      mn=input('Enter mode number: '); 
      if mn <1 | mn>nmodes
        disp('That mode doesn''t exist.  Try again.')
        mn=mnold;
      end
    end
    if answer=='b'
      clc
      scale=input('Enter scale factor: '); 
    end
	if answer=='d'
	  print
	end
  end

  if answer=='e'
    break
  end
  x=xx(:,mn);


%  y=zeros(ndof,1);
%  y(p)=x;x=y;
  for i=1:length(x)/3
    for j=1:2
      d(i,j)=node(i,j)+scale*x(i*3-3+j);
    end
  end
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
  hold off
  clf
  plot(node(:,1),node(:,2),'*b')
  uicontrol('style','pushbutton','units','normal','position',[.91 .95 .075 .05],'string','Print','callback','print')
  uicontrol('style','pushbutton','units','normal','position',[.91 .89 .075 .05],'string','Close','callback','delete(gcf)')
%%%
  dx=.01*(xh-xl);
  dy=.035*(yh-yl);
%  text(node(:,1)'+dx,node(:,2)'+dy,nnum)
  hold on
  plot(xh,yh,'w')
  plot(xl,yl,'w')
  plot(d(:,1),d(:,2),'or')
  text(d(:,1)'+dx,d(:,2)'+dy,nnum)
  ncons=size(ncon);
  for i=1:ncons(1)
    plot([node(ncon(i,1),1) node(ncon(i,2),1)],[node(ncon(i,1),2) node(ncon(i,2),2)],'--b')
  end
  for i=1:ncons(1)
    plot([d(ncon(i,1),1) d(ncon(i,2),1)],[d(ncon(i,1),2) d(ncon(i,2),2)],'r')
  end
  tinm1=['Mode ', num2str(mn)];
  tinm2=[', ' num2str(f(mn)/2/pi) , 'Hz'];
  titlename=[tinm1 tinm2];
  title(titlename)
  axis('image')
end
grid on

%Automatically check for updates
vtbchk

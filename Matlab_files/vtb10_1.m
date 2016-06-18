function vtb10_1(a,b,dt,len,n,func)
%VTB10_1 VTB10_1(xr,dxr,dt,len,n,func) plots the phase plane of the 
%    second order single variable differential equation given 
%    by func.  The variable 'func' is the name of the 'm' file
%    which gives yields the second time derivative of the variable
%    'x' given 'x' and the first time derivative of 'x'.  For 
%    example, see the example file vdp.m which is used to plot 
%    the phase plane of the Van der Pol equation.  The equation
%    must be in element by element form (i.e. './' instead of '/')
%    since it will be sent multiple values of position and 
%    velocity simultaneously in matrix form.
%    INPUT ARGUMENTS:
%    xr  : The range and step size of starting points for x.
%          i.e. [-3:.1:3] to start from x values from -3 to 3 by
%          step sizes of .1.
%    dxr : The range and step size of starting points for the 
%          velocity.
%    dt  : Step size of time for calculation of time derivatives.
%          Recommended step size is on the order of 1e-6 times
%          the total range for x.  The smaller, the more accurate
%          the plot.
%    len : Length of line drawn at each time step.  The smaller,
%          the more accurate.  However, 1e-2 time the total range
%          for x should give quite good results for simple phase
%          plots.
%    n   : Number of time steps to be taken.  More time steps
%          must be taken to truly see the phase portrait.
%          More also take more time linearly.
%    func: The name of the function which yields the second time
%          derivative of x given x and the first time derivative
%          of x.  It must accept only the inputs x and the
%          derivative of x and in that order.  See the function
%          'VTB10_ex.m'
%
%    Example:  Enter the command
%    vtb10_1([-3:1.5:3],[-3:1.5:3],.0001,.04,100,'vtb10_ex')
%
%    Pendulum Example:
%    vtb10_1([-3:1.5:3],[-3:1.5:3],.001,.04,200,'vtb10_expend')

%          Copyright Joseph C. Slater, April 19, 1992
%          
%          
%
nzero=0;
clf
[x,xd]=meshgrid(a,b);
s=size(x);
axis('square')
%axis([x(1,1) x(1,s(2)) xd(s(1),1) xd(1,1)]);

aa=version;ll=length(aa);

plot(x(1,1),xd(1,1))
hold on
range=max(a)-min(a);
for k=1:n
  x2=xd*dt+x;
  xdd=feval(func,x,xd);
  xd2=xd+xdd*dt;
  l=((x2-x).^2+(xd2-xd).^2).^.5;
  x2=(x2-x)./(l+1e-10*range)*len+x;
  xd2=(xd2-xd)./(l+1e-10*range)*len+xd;
  for i=1:s(1)
     for j=1:s(2)
        if xd(i,j)==0 & xdd(i,j)==0
           if exist('zlog')==0
              xx=num2str(x(i,j));
              xxd=num2str(xd(i,j));
              text1=['The point (' xx ',' xxd ') has a zero'];
              text2=[' first and second derivative.'];
              text1=[text1 text2]; 
              plot(x(i,j),xd(i,j),'o')
              disp(text1)
              nzero=nzero+1;
              zlog(nzero,1:2)=[x(i,j) xd(i,j)];
			  plot(x(i,j),xd(i,j),'or')
             else		  
              %[ones(l,1)*b(1) ones(l,1)*b(2)];
              old=sum((zlog'==[ones(nzero,1)*x(i,j) & ones(nzero,1)*xd(i,j)]'));
              if ~old
                 xx=num2str(x(i,j));
                 xxd=num2str(xd(i,j));
                 text1=['The point (' xx ',' xxd ') has a zero'];
                 text2=[' first and second derivative.'];
                 text1=[text1 text2]; 
                 disp(text1)
                 nzero=nzero+1;
                 zlog(nzero,1:2)=[x(i,j) xd(i,j)];
				 plot(x(i,j),xd(i,j),'or')
              end			  
           end
        end
        plot([x(i,j) x2(i,j)],[xd(i,j) xd2(i,j)])
     end
  end
  x=x2;
  xd=xd2;
  %drawnow
end
plot(zlog(:,1),zlog(:,2),'or')
hold off
grid on

%Automatically check for updates
vtbchk

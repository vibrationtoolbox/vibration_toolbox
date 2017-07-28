function vtb1_1(m,c,k,x0,v0,tf)

%VTB1_1 Free response of a single degree of freedom system.
% VTB1_1(m,c,k,x0,v0,tf) plots the free response of a single degree 
% of freedom system.  The arguments x0 and v0 represent the initial 
% conditions and tf represents the total time of the response. The
% system is described by mass m, damping c, and stiffness k.
% VTB1_1(zeta,w,x0,v0,tf) plots the free response of a single degree 
% of freedom system.  The arguments x0 and v0 represent the initial 
% conditions and tf represents the total time of the response.
% The system is in non-dimensional form, where zeta is the damping
% ratio, and w is the natural frequency in rad/s.
% Example: Free response for m=1, c=.01, k=2, x0=1, v0=0. For 100 sec.
% Type the following on the command line after the prompt:
% vtb1_1(1,.01,2,1,0,100)

%

figure
%uicontrol('style','pushbutton','units','normal','position',[.91 .95 .075...
% .05],'string','Print','callback','print')
%uicontrol('style','pushbutton','units','normal','position',[.91 .89 .075...
% .05],'string','Close','callback','delete(gcf)')

% This loop determines which type of input format you are using.
if nargin==5
  z=m;w=c;tf=v0;v0=x0;x0=k;m=1;c=2*z*w;k=w^2;
end

w=sqrt(k/m);
z=c/2/w/m;%(1.30)
wd=w*sqrt(1-z^2);%(1.37)

fprintf('The natural frequency is %.3g rad/s.\n',w);
fprintf('The damping ratio is %.3g.\n',z);
fprintf('The damped natural frequency is %.3g.\n',wd);
t=0:tf/1000:tf;
if z < 1
    A=sqrt(((v0+z*w*x0)^2+(x0*wd)^2)/wd^2);%(1.38)
    phi=atan2(x0*wd,v0+z*w*x0);%(1.38)
    x=A*exp(-z*w*t).*sin(wd*t+phi);%(1.36)
    fprintf('A= %.3g\n',A);
    fprintf('phi= %.3g\n',phi);
  elseif z==1
    a1=x0;%(1.46)
    a2=v0+w*x0;%(1.46)
    fprintf('a1= %.3g\n',a1);
    fprintf('a2= %.3g\n',a2);
    x=(a1+a2*t).*exp(-w*t);%(1.45)
  else
    a1=(-v0+(-z+sqrt(z^2-1))*w*x0)/2/w/sqrt(z^2-1);%(1.42)
    a2=(v0+(z+sqrt(z^2-1))*w*x0)/2/w/sqrt(z^2-1);%(1.43)
    fprintf('a1= %.3g\n',a1);
    fprintf('a2= %.3g\n',a2);
    x=exp(-z*w*t).*(a1*exp(-w*sqrt(z^2-1)*t)+a2*exp(w*sqrt(z^2-1)*t));%(1.41)
end
plot(t,x)
xlabel('Time')
ylabel('Displacement')
title('Displacement versus Time')

grid on

zoom on


%Automatically check for updates
vtbchk

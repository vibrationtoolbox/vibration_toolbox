function [t,x,v] = vtb1_3(rkf,u,t,x0,v0)
%VTB1_3   VTB1_3(xdd,f,t,x0,v0)
%       Runge-Kutta fourth order solution to a first order DE
%       t is a row vector from the initial time to the final time
%       by some time step. This time step should be approximately
%       1/10th the period of the system or less.  
%       x0 is the initial displacement.
%       v0 is the initial velocity.
%       The force vector f is ordered such that the nth column or row
%       of f is the force vector f evaluated at time n*dt.
%       xdd is a string variable of the the differential equation solved 
%           for the acceleration using x, v, t and f for the displacement,
%           velocity, time, and force. 
%       Example 1:
%          If the equation is
%          10 xdd + 2 xd + 1000 x= 30 sin(10 t)
%          then
%          xdd='-100*x-.2*v+3*sin(10*t)'
%
%       If the forcing is defined by a vector of numbers then 
%       Example 2:
%          If the equation is
%          10 xdd + 2 xd + 1000 x-x^3= 30 f
%          then
%          xdd='-100*x+0.1*x^3-.2*v+3*f'
%
%       Example 3 (Example 2 implemented):
%          t=0:.1:20;  % Creates time vector
%          f=0*t;      % Must be defined, even if zero
%          x0=1;  % Creates initial displacement
%          v0=0;  % Initial velocity
%          xdd='-100*x+0.1*x^3-.2*v+3*f';
%          [t,x,v]=vtb1_3(xdd,f,t,x0,v0); % Runs analysis.
%          plot(t,x); % Plots displacement versus time.
%          plot(t,v); % Plots velocity versus time.
%
%       Example 4:
%          t=0:.1:20;  % Creates time vector
%          f=0*t;      % Must be defined, even if zero
%          x0=1;  % Creates initial displacement
%          v0=0;  % Initial velocity
%          [t,x,v]=vtb1_3('-x-v/10+f/10+sin(2*t)',f,t,x0,v0); % Runs analysis.
%          plot(t,x); % Plots displacement versus time.
%          plot(t,v); % Plots velocity versus time.
%
%       Executing vtb1_3 without input arguments runs example 3.
%
%  This code WILL NOT give as good results as Matlab's integrators. It is
%  designed to be very simple to execute. It attempts to adjust the time
%  step for best results, but the user should try smaller time steps until 
%  convergence to make sure the results are reasonable. 
%

if nargin==0
   disp('Running demo 3')
   t=0:.1:20;  % Creates time vector
   f=0*t;      % Must be defined, even if zero
   x0=1;  % Creates initial displacement
   v0=0;  % Initial
   [x,v]=vtb1_3('-x-v/10+f/10',f,t,x0,v0); % Runs analysis.
   plot(t,x); % Plots displacement versus time.
   xlabel('Time')
   ylabel('Displacement')
   title('Displacement versus time, vtb1\_3 Example 3')
   grid on
   pause
   plot(t,v); % Plots velocity versus time.
   xlabel('Time')
   title('Velocity versus time, vtb1\_3 Example 3')
   ylabel('Velocity')
   grid on
   clear x v
else

    
xis=min(strfind(rkf,'*x'));
s=rkf(1:(xis-1));
lockom=min([strfind(s,'+') strfind(s,'-')]);
rkf((lockom+1):xis);
kom=str2num(rkf((lockom+1):(xis-1)));
if isempty(kom)
    kom=1;
end

TT=2*pi/sqrt(kom);
h=t(2)-t(1);

if h>TT/20
    disp(['You should probably have uses a smaller time step. It has been reduced to ' num2str(TT/20) '.'])
    x=[];v=[];
    if norm(u)>1e-14
    return
    else
    t=min(t):(TT/20):max(t);
    u=t*0;
    end
end

n=length(t);
%if n<100;
%    t=min(t):(max(t)-min(t))/500:max(t);
%end

z=zeros(2,length(t));
z(:,1)=[x0;v0];
h=t(2)-t(1);
tt=t;
for l1=1:(n-1);
   z1=z(:,l1);
   xx=z1(1);vv=z1(2);
   u1=u(l1);
   u2=u(l1+1);
   t=tt(l1);
   f=u1;
   x=xx;
   v=vv;
   k1=h*[v;eval(rkf)];
%   k1=h*feval(rkf,z1,u1,t(l1))

   x=xx+k1(1)/2;
   v=vv+k1(2)/2;
   k2=h*[v;eval(rkf)];
%   k2=h*feval(rkf,z1+.5*k1,u1,t(l1)+.5*h)

   x=xx+k2(1)/2;
   v=vv+k2(2)/2;
   k3=h*[v;eval(rkf)];
%   k3=h*feval(rkf,z1+.5*k2,u1,t(l1)+.5*h)

   x=xx+k3(1);
   v=vv+k3(2);
   k4=h*[v;eval(rkf)];
%   k4=h*feval(rkf,z1+k3,u1,t(l1)+h)
      %pause

%    k1=h*feval(rkf,z1,u1,t(l1));
%    k2=h*feval(rkf,z1+.5*k1,u1,t(l1)+.5*h);
%    k3=h*feval(rkf,z1+.5*k2,u1,t(l1)+.5*h);
%    k4=h*feval(rkf,z1+k3,u1,t(l1)+h);
   z(:,l1+1)=z(:,l1)+1/6*(k1+2*k2+2*k3+k4);
end

x=z(1,:);
v=z(2,:);
t=tt;
end

%Automatically check for updates
vtbchk

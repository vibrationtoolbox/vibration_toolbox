function [x,xd]=VTB9_2(n,dt,x0,xd0,a,b,c)

%VTB9_2 Finds the multiple degree of freedom system
%  undamped response using Euler's method.  
%  x=VTB9_2(n,dt,x0,a)  Solves the system given the 
%    initial state vector 'x0', the state matrix 'a', 
%    the time step to be used 'dt', and the number of 
%    steps to take 'n'.
%  [x,v]=VTB9_2(n,dt,x0,v0,m,d,k)  Solves the system 
%    given the initial displacement vector 'x0',
%    the mass matrix 'm', the stiffness matrix 'k', and 
%    the damping matrix 'd'. x0 is the initial displacement 
%    and v0 is the initial velocity. The remaining parameters 
%    are as described above.
%  The outputs are in the form of a matrix where each column 
%  represents the states a one time step and the rows represent 
%  a state as a function of time.
%       Example
%       a=[0 1;-1 -.1]; % Creates state matrix.
%       dt=.05;  % Time step size.
%       n=80;  % Number of steps to take.
%       x0=[1 ;0];  % Creates initial state vector.
%       [x,v]=vtb9_2(n,dt,x0,a); % Runs analysis.
%       t=0:dt:dt*n; % Creates time vector.
%       plot(t,x); % Plots displacement verses time.
%       pause
%       plot(t,v); % Plots velocity verses time.

disp('VTB9_2 has been grandfathered. Please use VTB1_2 in the future.')


if nargin==7
  la=length(a);
  x0=[x0;xd0];
  a=[zeros(la) eye(la);-a\c -a\b];
 else
  a=xd0;
end

% -----------------------------------------------
%             Where the action be at.
% -----------------------------------------------
x(:,1)=x0;
for i=2:1:n+1
  x(:,i)=a*dt*x(:,i-1)+x(:,i-1);
end
% -----------------------------------------------
%             Where the action be finished.
% -----------------------------------------------

lx=length(a);
if lx/2==floor(lx/2)
  xd=x(lx/2+1:lx,:);
  x=x(1:lx/2,:);
 else
  xd=x;
end

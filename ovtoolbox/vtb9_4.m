function [x,xd]=VTB9_4(n,dt,x0,xd0,a,b,c,u)

%VTB9_4 Solves the forced multiples degree of freedom system using Euler's 
%  method.  
%  x=VTB9_4(n,dt,x0,a,u)  Solves the system given the initial state vector 'x0',
%    the state matrix 'a', the time step to be used 'dt', and the number of 
%    steps to take 'n'.  The force matrix u is ordered such that the nth
%    column of u is the force vector u evaluated at time n*dt.
%  [x,xd]=VTB9_4(n,dt,x0,v0,m,d,k,u)  Solves the system given the initial 
%    displacement vector 'x0',
%    the mass matrix 'm', the stiffness matrix 'k', and the damping matrix 'd'.
%    The remaining parameters are as described above.
%  The outputs are in the form of a matrix where each column represents the 
%  states a one time step and the rows represent a state as a function of time.

% 	J. C. Slater 1-10-93
%	Revised J. C. Slater 8-17-94
%	Copyright (c) 1993-94 by Joseph C. Slater

disp('VTB9_4 has been grandfathered. Please use VTB1_4 in the future.')


if nargin==8
   la=length(a);
   x0=[x0;xd0];
   a=[zeros(la) eye(la);-a\c -a\b];
   u=[u*0;u];
  else
   u=a;
   a=xd0;
end
 

% -----------------------------------------------
%             Where the action is.
% -----------------------------------------------
x(:,1)=x0;
for i=2:1:n+1
  x(:,i)=a*dt*x(:,i-1)+x(:,i-1)+dt*u(:,i-1);
end
% -----------------------------------------------
%             Where the action is finished.
% -----------------------------------------------

lx=length(a);
if lx/2==floor(lx/2)
  xd=x(lx/2+1:lx,:);
  x=x(1:lx/2,:);
 else
  xd=x;
end

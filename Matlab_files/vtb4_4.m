function [x,xd]=vtb4_4(M,D,K,F,T, V0, X0)

%VTB1_4 Solves the forced multiple degree of freedom system matrix
%exponential methods (closed-form accurate for linear systems)
%using Matlab's lsim.
%  
%  [x,xd]=vtb4_4(M,C,K,F,T, V0, X0)  Solves the system given the initial 
%    displacement vector 'X0', initial velocity vector 'V0',
%    the mass matrix 'M', the stiffness matrix 'M', and the damping
%    matrix 'C' and force 'F'.
%    T is a row vector of evenly spaced times. 
%    F is a matrix of forces over time, each column corresponding
%    to the corresponding column of T, each row corresponding to
%    the same numbered DOF.
%  
%  Time step should be chosen to be much smaller than the smallest
%  period of your system. 
%
%  EXAMPLE (4.6.1, 3rd edition)
%  M=[9 0;0 1]
%  K=[27 -3;-3 3]
%  C=K/10
%  x0=[0;1] %Making these up (not given)
%  v0=[1;0] %Again making these up
%  t=0:.1:100; % The number in the middle should be 1/10 the
%              %smallest period
%  F=[0*t;3*cos(2*t)];% Don't forget the zero forces on row 1
%  [x,xd]=vtb4_4(M,C,K,F,t,v0,x0);  
%  plot(t,x(1,:),t,.2451*cos(2*t-2.9442)-.6249*sin(2*t)) %plots
%  %response of DOF 1, compares to closed form solution
% 
  
% 	J. C. Slater 5/29/09
%	Copyright (c) 2009 by Joseph C. Slater

A=[zeros(size(M)), eye(size(M));-M\K,-M\D];

B=[zeros(size(M));inv(M)];

C=[eye(2*size(M))];

D=0*B;

F=F';
T=T';
sys=ss(A,B,C,D);

y=lsim(sys,F,T,[X0;V0],'foh');



x=y(:,1:size(M,1))';
xd=y(:,[1:size(M,1)]+size(M,1))';

  


%Automatically check for updates
vtbchk

function [v,w,zeta]=VTB9_1(m,d,k)
%VTB9_1 VTB9_1(m,d,k)
% [v,w,zeta]=VTB9_1(m,d,k)
% function VTB9_1 finds the mode shapes and natural frequencies of
% a linear second order matrix equation.  
% [v,w]=VTB9_1(m,k) finds the mode shapes and natural frequencies 
% for the undamped case.

disp('VTB9_1 has been grandfathered. Please use VTB4_3 in the future.')

if nargin==2
  k=d;
  [v,w]=eig(m\k);
  w=sqrt(w);
end
if nargin==3
  if norm(d/m*k-k/m*d) < 1e-8*norm(k/m*d)
    disp('Damping is proportional, eigenvectors are real.')
    [v,w]=eig(m\k);
    w=sqrt(w);
    zeta=(v'*m*v)\(v'*d*v)/2/w;
   else
    disp('Damping is non-proportional, eigenvectors are complex.')
    a=[0*k eye(length(k));-m\k -m\d];
    [v,w1]=eig(a);
    w=abs(w1);
    zeta=-real(w1)/w;
  end
end
w=diag(w);zeta=diag(zeta);

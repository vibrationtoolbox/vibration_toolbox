function [v,w,zeta]=vtb4_3(m,d,k)
%VTB4_3 VTB4_3(m,d,k)
% [u,w,zeta]=VTB4_3(m,d,k)
% function VTB4_3 finds the mode shapes and natural frequencies of
% a linear second order matrix equation.  
% [u,w]=VTB4_3(m,k) finds the mode shapes and natural frequencies 
% for the undamped case.
if nargin==2
  k=d;
  [v,w]=eig(m\k);
  w=sqrt(diag(w));
  [d,i]=sort(sqrt(w/2/pi));
  w=w(i);
  v=real(v(:,i));
  vnorms=sqrt(v'*m*v);	
  v=real(v/vnorms);
end
if nargin==3
  if norm(d/m*k-k/m*d) < 1e-8*norm(k/m*d)
    disp('Damping is proportional, eigenvectors are real.')
    [v,w]=eig(m\k);
    w=sqrt(diag(w));
    [f,i]=sort(sqrt(w/2/pi));
    w=w(i);
    v=real(v(:,i));
    vnorms=sqrt(v'*m*v);	
    v=real(v/vnorms);
    zeta=diag((v'*m*v)\(v'*d*v)/2/diag(w));
   else
    disp('Damping is non-proportional, eigenvectors are complex.')
    a=[0*k eye(length(k));-m\k -m\d];
    [v,w1]=eig(a);
    w=abs(w1);
    zeta=-real(w1)/w;
  end
end

%if nargout==3;zeta=diag(zeta);end

%Automatically check for updates
vtbchk

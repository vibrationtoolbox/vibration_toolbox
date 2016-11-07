function [lambda,phi]=vtb7_6(M,K,n)
%  VTB7_6(M,K)
%
%  [lambda,phi] = VTB7_6(M,K) produces returns the first eigenvalue
%  lambda=(omega^2) and mode shape, phi using the power method. 

% Copyright Joseph C. Slater, 7/14/2007
% Shift eigenvalues to allow solution for rigid body modes. 

if nargin==2
	n=1;
end
shift=1;
K=K+shift*M;
delta=10;
X=randn(size(M,1),1);
X=X/norm(X);
%X=[0;1]
oldev=0;
A=K\M;
i=0;

for i=1:n

	while delta>100*eps
		Xnew=A*X;
		newev=norm(Xnew);
		delta=abs(newev-oldev);
		X=Xnew/norm(Xnew);
		oldev=newev;
	end
	lambda(i)=1/newev-shift;
	phi(:,i)=X;
	A=matrdefl(A,newev,X);
	delta=10;
end

%phi=X;
%lambda=1/newev-shift;



function Adef=matrdefl(A,lam,phi)
Adef=A-lam*phi*phi';

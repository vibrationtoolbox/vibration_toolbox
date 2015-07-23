function [x,v,t]=vtb4_2(M,K,x0,v0,tf,plotpar)
%VTB4_2 Free response of an undamped system.
% [x,v,t]=VTB4_2(M,K,x0,v0,tf) will return the displacement (x) and
% velocity (v) response for the system described by the mass
% matrix M and stiffness matrix K.  The vector of time points, t,
% is also returned.  The initial displacement, x0, and initial 
% velocity, v0, must also be a column vector input into the 
% function. The response is plotted from zero to final time tf.
% [x,v,t]=VTB4_2(M,K,x0,v0,tf,plotpar) will also plot the various 
% responses depending on the value of plotpar:
%
% plotpar = 1: plot only displacements (x)
% plotpar = 2: plot only velocities (v)
% plotpar = 3: plot both displacements and velocities

if nargin == 5
  plotpar=0;
end
% if plotpar~=0
% 
%   aa=version;ll=length(aa);
%   
% 
% end
t=0:tf/200:tf;

N=length(M);

A=[zeros(N) eye(N);-M\K zeros(N)];
%B=zeros(2*N,1);

IC=[x0;v0];

Ad=expm(A*(t(2)-t(1)));
eig(Ad);

X(:,1)=IC;

for i=2:length(t)
   X(:,i)=Ad*X(:,i-1);
end

x=X(1:N,:);
v=X(N+1:2*N,:);

%plotting routines

if plotpar~= 0
   if plotpar~=2
      for i=1:N
        plot(t,x(i,:))
        ylabel(['Displacement of X',num2str(i)])
        xlabel('time (sec)')
        title('Press any key to continue')
        grid on
		pause
      end
   end 
   if plotpar~=1
      for i=1:N
        plot(t,v(i,:))
        ylabel(['Velocity of X',num2str(i)])
        xlabel('time (sec)')
        title('Press any key to continue')
        grid on
		pause
      end
   end
end


%Automatically check for updates
vtbchk

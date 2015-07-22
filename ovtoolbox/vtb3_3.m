function [ap,bp]=vtb3_3(dat,t,n)
%VTB3_3  Fourier series approximation to a function.
%[a,b]=VTB3_3(dat,t,n) returns Fourier coefficients of a function
%  The coefficients are numerical approximations of the true
%  coefficients.
%    dat is a vector of data representing the function
%    t   is the corresponding time vector
%    n   is the desired number of terms to use in the Fourier series
% Intermediary plots show the impact of each successive term 
% on the total seres if no output arguments are specified.
% EXAMPLE: Manually performs the steps of the command vtb3_3(5)
% f=[ -1:.04:.96 1:-.04:-.96]'+1;
% t=(0:length(f)-1)/length(f)';
% plot(t,f)
% [a,b]=vtb3_3(f,t,5);
% vtb3_3(f,t,5)
%
% VTB3_3(N) displays the N term Fourier approximation to a 
% triangular input.  The approximation is plotted versus time 
% normalized by the period of the wave. 
%
% VTB3_3 displays the 5 term Fourier approximation to a 
% triangular input.  The approximation is plotted versus time 
% normalized by the period of the wave. 
%
% Note that these results are only an approximation, and the quality
% depends on the number of points used, and the proper selection of begining
% and end points. 


% Copyright Joseph C. Slater, Dec 1996
% Revised 03/13/03 - Now runs on Octave.
% Revised 02/29/00 - Now can run with no arguments.
% Revised 11/11/98 - Example changed to match default function 
%                    (Example 3.3.1)
% Revised 12/10/97 - Improved location of legend to avoid covering up data.
%                    Disclaimer on quality
clg
clc
if nargin==0
  if nargout==0
	vtb3_3(5);
  else
        [a,b]=vtb3_3(5);
  end
else

if nargin==1
  n=dat;
  tau1=0:.01:.5;
  Ftr1=(4*tau1-1);
  tau2=.51:.01:.99;
  Ftr2=3-4*tau2;
  t=[tau1 tau2]';
  dat=[Ftr1 Ftr2]';
end

if size(dat,1)==1
  dat=dat';
end

if size(t,1)==1
  t=t';
end

len=length(dat)/2;

grid('on')
fs=(fft(dat))/len;
fs(1:10);
a0=fs(1);

a=[a0; real(fs(2:length(fs/2)))];
b=-imag(fs(2:length(fs/2)));
len=len*2;
dt=2*pi/len;
tp=(0:dt:2*pi-dt)';

datapprox=a(1)/2+zeros(size(dat));

plot(t,dat,"-;Data;",t,datapprox,'o;New Approximation;')


if nargout==0
  context=['Press Return to continue. i = ' num2str(i) '.'];
  disp(context)
  pause
end

for i=2:n+1
  hold off

  newdat=a(i)*cos(tp*(i-1))+b(i-1)*sin(tp*(i-1));
  datapprox=datapprox+newdat;
  if nargout==0
    plot(t,dat,"-;Data;",t,datapprox,...
	 ['o;' num2str(i-2) ' Term Approximation;'],t,datapprox-newdat,...
	 ['x;' num2str(i-1) ' Term Approximation;'],t,newdat,...
	 ["+;Contribution of New Terms (i = " num2str(i-1) ");"])
    context=['Contribution of terms n=' num2str(i-1)];
    disp(['Press Return to continue. i = ' num2str(i-1) '.'])
    pause
  end
end
if nargout==2
ap=a(1:n+1);
bp=b(1:n);
end
%a(1:3)
%b(1:3)

nargout;
%if nargout~=0
%  ap=a(1:n+1);bp=b(1:n);
%end
%end
%a
%b
end
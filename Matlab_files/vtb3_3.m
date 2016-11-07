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
% EXAMPLE: Find the Fourier series of a square wave. It nicely demonstrates
% the Gibb's effect. 
% t=0:.001:1;
% f=2*(t<0.5);
% [a,b]=vtb3_3(f,t,15);
% vtb3_3(f,t,15)
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
% Revised 10/20/11 - Added example. 
% Revised 02/29/00 - Now can run with no arguments.
% Revised 11/11/98 - Example changed to match default function 
%                    (Example 3.3.1)
% Revised 12/10/97 - Improved location of legend to avoid covering up data.
%                    Disclaimer on quality

if nargin==0
  vtb3_3(5)
  % If this is a demo mode, we don't want to do anythin after
  % running vtb3_3(5)
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
  if nargin==2
    n=100;
    disp('Using default value of n=100')
  end
  if size(dat,1)==1
    dat=dat';
  end
  
  if size(t,1)==1
    t=t';
  end

  len=length(dat)/2;
  plot(dat),grid on
  fs=(fft(dat))/len;
  
  a0=fs(1);

  a=[a0; real(fs(2:length(fs/2)))];
  b=-imag(fs(2:length(fs/2)));
  len=len*2;
  dt=2*pi/len;
  tp=(0:dt:2*pi-dt)';
  
  datapprox=a(1)/2+zeros(size(dat));
  plot(t,dat,t,datapprox)
  grid on
  %aa=version;ll=length(aa);
  grid on
  %context=['Contribution of terms n=' num2str(i-1)];
  legend('Function','New Approximation')
  if nargout==0
    disp('Press ''enter'' to advance')
    pause
  end
  
  for ii=2:n+1
    %  a(ii)
    %  b(ii-1)
    newdat=a(ii)*cos(tp*(ii-1))+b(ii-1)*sin(tp*(ii-1));
    datapprox=datapprox+newdat;
    if nargout==0
      %legend off
      plot(t,dat,t,datapprox,'o',t,datapprox-newdat,'x',t,newdat)
      %pause
      %aa=version;ll=length(aa);
      grid on
      context=['Contribution of terms n=' num2str(ii-1)];
      legend('Function','New Approximation','Old Approximation',context,0)
      pause
    end
  end
  
  legend off
  plot(t,dat,t,datapprox),grid on
  
  legend('Function','Approximation')
  %aa=version;ll=length(aa);
  
  if nargout~=0
    ap=a(1:n+1);bp=b(1:n);
  end
end
%Automatically check for updates
vtbchk

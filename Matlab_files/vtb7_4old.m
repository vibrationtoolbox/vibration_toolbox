function [z,nf,a,com]=VTB7_4(f,TF,b)
%[z,nf,a,com]=VTB7_4(f,TF) Curve fit to SDOF FRF.
% f is the frequency vector in Hz. It does not have to 
%    start at 0 Hz.
% TF is the complex transfer function.
% z and nf are the damping ratio and natural frequency (Hz)
% a is the product of the residues of the coordinates the 
% transfer function is between. (For example, in the example 
% below, a1 times a2 is returned. See equation 7.42)
% If com is returned as a real number, then it is the 
% compliance between the two coordinates. 
% Only one peak may exist in the segment of the FRF passed to 
% VTB7_4. No zeros may exist withing this segment. Otherwise, 
% curve fitting becomes unreliable.
%
% EXAMPLE:
% M=eye(2);
% K=[2 -1;-1 2];
% C=.01*K;
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,2,linspace(0,.5,1024));
% figure(1)
% n=250;
% f2=Freq((1:n)+450);
% R2=Recep((1:n)+450);
% R2=R2+.1*randn(n,1)+.1*randn(n,1)*i;% Poorly Simulated Noise
% [z,nf,a,com]=vtb7_4(f2,R2)
%
% Note that by changing the parts of Freq and Recep used
% We can curve fit to other modes.

% Copyright Joseph C. Slater, 10/8/99
% Updated 11/8/99 to improve robustness


global XoF

if nargin==2
	[y,in]=max(abs(TF));
	lf=length(f);
	f(in);
	a0=abs(TF(1))*(2*pi*f(in))^2;
	z=.0005;
	a0=-sign(imag(TF(in)))*abs(TF(in))*2*z*(2*pi*f(in))^2;TF(in);
	x=[a0;z;2*pi*f(in);0;0;0];%sign(real(TF(1)))*
	%x2=x;%
	%cost=vtb7_4(x,f,TF)
	if in-3<1|in+2>length(f)
		disp('The peak response must be near the middle of the data')
		disp('Please center your peak and try again')
		break
	end
		
	x	
	x=fmins('vtb7_4',x,[],[],f(in-3:in+2),TF(in-3:in+2));
    x
	cferr(x,f,TF);	%x
	%cost=vtb7_4(x,f,TF)
	x=fmins('vtb7_4',x,[],[],f,TF);
	x
	%cost=vtb7_4(x,f,TF)
	x=fmins('vtb7_4',x,[],[],f,TF);
	x
	%cost=vtb7_4(x,f,TF)
	x=fmins('vtb7_4',x,[],[],f,TF);
	x
	%cost=vtb7_4(x,f,TF)
	x=fmins('vtb7_4',x,[],[],f,TF);
	x
	%cost=vtb7_4(x,f,TF)
	%x=x2
	z=x(2);om=x(3);
	%z,om
	if f(1)==0
		k=x(4)+x(1)/om^2;
	else
		k=sqrt(-1);
	end
	com=k;
	nf=om/2/pi;%*2*pi;
	a=x(1);
	%plot(f,20*log10(abs(XoF)),'g',f,20*log10(abs(TF)))
	%grid on
	%zoom on
	if 1==1
	  Fmin=min(f);
	  Fmax=max(f);
	  phase=unwrap(angle(TF))*180/pi;
	  phase2=unwrap(angle(XoF))*180/pi;size(phase);
		%size(XoF)
	  subplot(2,1,1)
	  plot(f,20*log10(abs(XoF)),f,20*log10(abs(TF)))
	  as=axis;
	  zoom on
	  legend('Identified FRF','Experimental FRF',0)
	  axis([Fmin Fmax as(3) as(4)])
	  xlabel('Frequency (Hz)')
	  ylabel('Mag (dB)')
	  grid on
	%  Fmin,Fmax,min(mag),max(mag)
	%  axis([Fmin Fmax minmag maxmag])
	  while phase2(in)>50
		  phase2=phase2-360;
	  end
	  phased=phase2(in)-phase(in);
	  phase=phase+round(phased/360)*360;
	  phmin_max=[floor(min(min([phase;phase2]))/45)*45 ceil(max(max([phase;phase2]))/45)*45];
	  subplot(2,1,2)
	  plot(f,phase2,f,phase)
	  xlabel('Frequency (Hz)')
	  ylabel('Phase (deg)')
	  legend('Identified FRF','Experimental FRF',0)
	
	  grid on
	  axis([Fmin Fmax  phmin_max(1) phmin_max(2)])
	  gridmin_max=round(phmin_max/90)*90;
	  set(gca,'YTick',gridmin_max(1):22.5:gridmin_max(2))
	  zoom on
	
  end
else
%	global  XoF
%f,TF,b
    x=f;
	f=TF;
	TF=b;
	w2=f*2*pi;
	lx=length(x);
	x(3)=abs(x(3));
	x(2)=abs(x(2));
	XoF=x(lx-2)+x(lx-1)*i*w2-x(lx)*w2.^2;
%	for j=1:(lx/3)-1
	XoF=XoF+x(1)./(-w2.^2+2*x(2)*w2*i*x(3)+x(3)^2);
%	end
	
	vtb74=norm(XoF-TF);
	z=vtb74;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cferr=cferr(x,f,TF)
global  XoF
w2=f*2*pi;
lx=length(x);

XoF=x(lx-2)+x(lx-1)*i*w2-x(lx)*w2.^2;
for j=1:(lx/3)-1
XoF=XoF+x(3*j-2)./(-w2.^2+2*x(3*j-1)*w2*i*x(3*j)+x(3*j)^2);
end

cferr=norm(XoF-TF);
%pause

%Automatically check for updates
vtbchk

function vtb5_3(mu,beta)
%VTB5_3 Normalized displacement of the primary mass for absorber design.
% VTB5_3(mu,beta) plots the normalized displacement of the primary mass
% for a vibration absorber design.  The variable mu is defined
% as the absorber mass divided by the primary mass (ma/m).  The 
% ratio of natural frequencies (wa/wp) is defined as beta.
 

mb=mu*beta^2;
b=beta;
%eqn=[1 0 -(4+2*mb) 0 (1+(2+mb)^2) 0 -(2+2*mb)];
c8=b^4;
c7=0;
c6=-2*b^2*(1+b^2+b^2*mu);
c5=0;
c4=b^2*(4+b^2+2*mu+2*b^2*mu+b^2*mu^2);
c3=0;
c2=-2*b^2*(1+mu);

eqn=[c8 c7 c6 c5 c4 c3 c2];
r=sort(roots(eqn));       %Solves for crossings of 1
r1=r(4);
r2=r(5);

wdr=(0:.01:r(6)*1.2)';          %normalized frequency wdr/wa

one=ones(length(wdr),1);

num=one-wdr.^2;

den1=(one+mu*beta^2*one-beta^2*wdr.^2);
den2=(one-wdr.^2);
den=den1.*den2-mu*beta^2*one;

f=num./den;

%A=(1-r.^2)./...
%((one+mu*beta^2*one-beta^2*r.^2).*(one-r.^2)-mu*beta^2*one);

j=[0 1];


aa=version;ll=length(aa);

plot(wdr,abs(f),wdr,one,[r1 r1],j,'--',[r2 r2],j,'--')
%axis([0 2 0 3]);
axis([0 r(6)*1.2 0 3]);
grid on
hold on

%Draws shaded lines
sh=r1:.1*(r2-r1):r2;
%for i=1:length(sh)
%   plot([sh(i) sh(i)],j,'--')
%end
tmp1='Normalized magnitude of primary mass for \mu = ';
title([tmp1,num2str(mu),' and \beta = ',num2str(beta)])
ylabel('|X_k/F_o|')
xlabel('normalized frequency \omega_{dr}/\omega_{a}')
text(.4,.8,'Useful operating range','Units','normalized')
text(.4,.75,[num2str(r1),' < \omega_{dr}/\omega_{a} < ',num2str(r2)],'Units','normalized')
set(1,'Units','default')

%[x,y]=meshdom(r1:.02:r2,0:.05:1);
[x,y]=meshgrid(r1:.02:r2,0:.05:1);
plot(x,y,'.m')
grid on
zoom on
hold off
%Automatically check for updates
vtbchk


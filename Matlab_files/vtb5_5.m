function vtb5_5(ma,m,ca,ka,k,rfin)
%VTB5_5 Normalized amplitude of the primary mass for a 
% damped vibration absorber design.
% VTB5_5(ma,m,ca,ka,k,rfin) plots the normalized amplitude
% for a system with the parameters:
%  ma    = absorber mass
%  m     = primary mass
%  ca    = absorber damping coefficient
%  k     = primary stiffness
%  ka    = absorber stiffness
%  rfin  = maximum normalized frequency for the plot

%calculate variables for equation 5.37
wa=sqrt(ka/ma);
wp=sqrt(k/m);
b=wa/wp;
mu=ma/m;
z=ca/(2*ma*wp);

r=0:.01:rfin;       %r=wdr/wp
r=r';

num=(2*z*r).^2+(r.^2-b^2).^2;

den1=((2*z*r).^2).*(r.^2-1+mu*r.^2).^2;
den2=(mu*r.^2*b^2-(r.^2-1).*(r.^2-b^2)).^2;

f=sqrt(num./(den1+den2));

aa=version;ll=length(aa);

plot(r,f)
%axis([0 rfin 0 max(f)*1.5]);
axis([0 rfin 0 2]);
grid on
zoom on

title('Normalized amplitude of the primary mass')
ylabel('|Xk/Fo|')
xlabel('normalized frequency - wdr/wp')
text(.2,.85,['Primary mass frequency = ',num2str(wp),' rad/s'],'units','normal')
text(.2,.8,['Frequency ratio (beta) = ',num2str(b)],'units','normal')
text(.2,.75,['Mass ratio (mu) = ',num2str(mu)],'units','normal')
text(.2,.70,['Absorber damping ratio = ',num2str(z)],'units','normal')

%Automatically check for updates
vtbchk

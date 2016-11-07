function vtb5_4(b)
%VTB5_4 Mass ratio versus normalized frequency for absorber design.
% VTB5_4(beta) plots the mass ratio as a function of the normalized
% frequency (w/wa) for the given value of beta.

mu=0:.01:1;
mu=mu';
one=ones(length(mu),1);

wl=(one+b^2*(one+mu))/(2*b^2);
wr=(sqrt(b^4*(one+mu).^2-2*b^2*(one-mu)+one))/(2*b^2);

w1=sqrt(wl+wr);
w2=sqrt(wl-wr);

axis([0 w1(101)+.2 0 1]);

aa=version;ll=length(aa);


plot(w1,mu,w2,mu,'-')
tmp1='Mass ratio versus system natural frequency for \beta = ';
title([tmp1,num2str(b)])
ylabel('mass ratio - \mu')
xlabel('normalized frequency \omega/\omega_a')
grid on
axis;

%Automatically check for updates
vtbchk

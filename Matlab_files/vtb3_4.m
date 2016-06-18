function vtb3_4(f)
%VTB3_4 Response spectrum for a SDOF system.
% VTB3_4(f) will display the response spectrum to a partial ramp 
% input(see Figure 3.13) for the system with natural frequency 
% f (in Hz). 

t=linspace(.001*4/f,10/f,200);
size(t)

[rt,ct]=size(t);
%Checks to make sure t is a column vector
if rt<ct
   t=t';
end

w=2*pi*f;

one=ones(length(t),1);

%Breaks 3.81 into two parts
Rs1=one./(w*t);
Rs2=sqrt(2*(1-cos(w*t)));

%Calculates 3.81
Rs=one+Rs1.*Rs2;



plot(t,Rs)
grid on
title(['Response spectrum of a SDOF system with f = ',num2str(f),' Hz']);
ylabel('Dimensionless maximum response - (xk/Fo)max')
xlabel('Rise time (t_1)')


%Automatically check for updates
vtbchk

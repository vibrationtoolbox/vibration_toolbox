function [w,x,U]=vtb6_3(n,bctype,bmpar,npoints)

%VTB6_3 Natural frequencies and mass normalized mode shape for an Euler-
% Bernoulli beam with a chosen boundary condition.
% [w,x,U]=VTB6_3(n,bctype,bmpar,npoints) will return the nth natural 
% frequency (w) and mode shape (U) of an Euler-Bernoulli beam.
% If n is a vector, return the coresponding mode shapes and natural
% frequencies.
% With no output arguments the modes are ploted.
% If only one mode is requested, and there are no output arguments, the
% mode shape is animated.
% The boundary condition is defined as follows:
%
% bctype = 1 free-free
% bctype = 2 clamped-free
% bctype = 3 clamped-pinned
% bctype = 4 clamped-sliding
% bctype = 5 clamped-clamped
% bctype = 6 pinned-pinned
%
% The beam parameters are input through the vector bmpar:
% bmpar = [E I rho A L];
% where the variable names are consistent with Section 6.5 of the 
% text.
%
%% Example: 20 cm long aluminum beam with h=1.5 cm, b=3 cm
%% Animate the 4th mode for free-free boundary conditions
% E=7.31e10;
% I=1/12*.03*.015^3;
% rho=2747;
% A=.015*.03;
% L=0.2;
% vtb6_3(4,1,[E I rho A L]);
%

% Copyright Joseph C. Slater, 2007
% Engineering Vibration Toolbox

if nargin==3
   npoints=100;
end

E=bmpar(1);
I=bmpar(2);
rho=bmpar(3);
A=bmpar(4);
L=bmpar(5);

len=[0:(1/(npoints-1)):1]';  %Normalized length of the beam

%Determine natural frequencies and mode shapes depending on the
%boundary condition.

if bctype==1
    desc=['Free-Free '];
    Bnllow=[0  0 4.73004074486 7.8532046241 10.995607838 14.1371654913 17.2787596574];
    for i=1:length(n)
        if n(i)>7
            %for i=6:n
            Bnl(i)=(2*n(i)-3)*pi/2;
            %end
        else
            Bnl(i)=Bnllow(n(i));

        end
    end
    for i=1:length(n)
        if n(i)==1
            w(i,1)=0;
            U(:,i)=1+len*0;
        elseif n(i)==2
            w(i,1)=0;
            U(:,i)=len-.5;
        else
            sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
            w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
            b=Bnl(i)*len;
            U(:,i)=cosh(b)+cos(b)-sig*(sinh(b)+sin(b));
        end
        
        %U(:,i)=U(:,i)/U(101,i);
    end
    

elseif bctype==2
    desc=['Clamped-Free '];
    Bnllow=[1.88 4.69 7.85 10.99 14.14];
    for i=1:length(n)
        if n(i)>5
            %for i=6:n
            Bnl(i)=(2*n(i)-1)*pi/2;
            %end
        else
            Bnl(i)=Bnllow(n(i));
        end
    end
    for i=1:length(n)
        sig=(sinh(Bnl(i))-sin(Bnl(i)))/(cosh(Bnl(i))-cos(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(101,i);
    end
    
elseif bctype==3
    desc=['Clamped-Pinned '];
    Bnllow=[3.93 7.07 10.21 13.35 16.49];
    for i=1:length(n)
        if n(i)>5
            %for i=6:n
            %Bnl(i)=(2*n(i)-1)*pi/2
            Bnl(i)=(4*n(i)+1)*pi/4;
            %end
        else
            Bnl(i)=Bnllow(n(i));
        end
    end

    for i=1:length(n)
        sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(52,i);
    end

elseif bctype==4
    desc=['Clamped-Sliding '];
    Bnllow=[2.37 5.50 8.64 11.78 14.92];

    for i=1:length(n)
        if n(i)>5
            %for i=6:n
            %Bnl(i)=(2*n(i)-1)*pi/2
            Bnl(i)=(4*n(i)-1)*pi/4;
        else
            Bnl(i)=Bnllow(n(i));
        end
    end
    for i=1:length(n)
        sig=(sinh(Bnl(i))+sin(Bnl(i)))/(cosh(Bnl(i))-cos(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(101,i);
    end

elseif bctype==5
    desc=['Clamped-Clamped']
    Bnllow=[4.73 7.85 11 14.14 17.28];
    
        for i=1:length(n)
        if n(i)>5
            %for i=6:n
            %Bnl(i)=(2*n(i)-1)*pi/2
            Bnl(i)=(2*n(i)+1)*pi/2;
        else
            Bnl(i)=Bnllow(n(i));
        end
    end

    for i=1:length(n)
        sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        b=Bnl(i)*len;
        U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
        %U(:,i)=U(:,i)/U(52,i);
    end
elseif bctype==6
    desc=['Pinned-Pinned'];
    for i=1:length(n)
        Bnl(i)=n(i)*pi;
        w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
        U(:,i)=sin(Bnl(i)*len);
    end
end

for i=1:length(n)
    U(:,i)=U(:,i)/sqrt(U(:,i)'*U(:,i)*rho*A*L);   
end


%stopstop=0;pausepause=0
global stopstop ppause;
ppause=0;
x=len*L;
%Plotting routine if so chosen.
if nargout==0
    if length(n)~=1
        for i=1:length(n)
            plot(x,U(:,i))
            axis([0 L min(min(U)) max(max(U))])
            figure(gcf)
            title([desc,'  ','Mode ',int2str(i),'     Natural Frequency = ',num2str(w(i)),' rad/s'])
            ylabel('Modal Amplitude')
            xlabel('Length along bar - x')
            grid on
            disp('Press return to continue')
            pause
        end
    else
        nsteps=50;
        clf
        step=2*pi/(nsteps);
        i=0:step:(2*pi-step);
        hold off
        handle=uicontrol('style','pushbutton','units','normal','backgroundcolor','red','position', ...
            [0.94 .94 .05 .05],'String','Stop','callback','global stopstop;stopstop=1;');
        handle2=uicontrol('style','pushbutton','units','normal','backgroundcolor','yellow','position', ...
            [0.94 .87 .05 .05],'String','Pause','callback','global ppause;ppause=1;');
        handle3=uicontrol('style','pushbutton','units','normal','backgroundcolor','green','position', ...
            [0.94 .80 .05 .05],'String','Resume','callback','global ppause;ppause=0;');
        
        stopstop=0;
        bb=0;
        while stopstop==0&bb<100
            bb=bb+1;
            for ii=[i  ]
                while ppause==1
                    pause(.01)
                    if stopstop==1
                        delete(handle), delete(handle2), delete(handle3)
                        return
                    end
                    
                end
                
                plot(x,U(:,1)*cos(ii))
                axis([0 L -max(abs(U)) max(abs(U))])
                grid on
                figure(gcf)
                title([desc,'  ','Mode ',int2str(n),'     \omega_n = ',num2str(w(1)),' rad/s'])
                ylabel('Modal Amplitude')
                xlabel('Length along bar - x')
                drawnow
                %pause
            end
        end
        clear stopstop
        delete(handle), delete(handle2), delete(handle3)
    end
end

%Automatically check for updates
vtbchk

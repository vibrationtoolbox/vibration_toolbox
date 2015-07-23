function vtb3_5(a0,aodd,aeven,bodd,beven,N,T,o)
%VTB3_5  Plot Fourier Series Approximation of Function.
%  VTB3_5(a0,a,b,N,T,o) Plot the Fourier series defined by:
%  N is the number of terms. Other arguments as defined by Inman's text. 
%  a and b should ideally be strings as shown in the example.
%  VTB3_5(a0,aodd,aeven,bodd,beven,N,T,o) Plot the Fourier with differing
%  even and odd values. 
%
%  Setting the optional variable o=1 disables value overrides. This can 
%  cause run times to be excessive in some cases. 
%
%  Example 1:
%  %Plot a 20 term Fourier series representation of a shifted square wave
%  %High value should be 1, low value should be -2.
%
%  vtb3_5(-1,0,'-3*(-1+(-1)^n)/n/pi',20,2)
%
%  Example 2:
%  Example 3.3.1 (3rd edition)
%
%  vtb3_5(0,'-8/pi^2/n^2',0,0,0,20,10)


% Copyright Joseph C. Slater, Dec 2006
% Revised Nov, 2007 to simplify use with differing even and odd values

if nargin>0
	if nargin<6||nargin==7
		o=0;
	end
	if nargin<7&&nargin>0
		if ~ischar(a0)
			a0=num2str(a0);
		end

		if ~ischar(aodd)
			a=num2str(aodd);
		else
			a=aodd;
		end

		if ~ischar(aeven)
			b=num2str(aeven);
		else
			b=aeven;
		end

		if ischar(bodd)
			N=str2double(bodd);
		else
			N=bodd;
		end

		if ischar(beven)
			T=str2double(beven);
		else
			T=beven;
		end
	elseif nargin>6
		if ~ischar(a0)
			a0=num2str(a0);
		end

		if ~ischar(aodd)
			aodd=num2str(aodd);
		end

		if ~ischar(aeven)
			aeven=num2str(aeven);
		end

		if ~ischar(bodd)
			bodd=num2str(bodd);
		end

		if ~ischar(beven)
			beven=num2str(beven);
		end
		
		if ischar(T)
			T=str2double(T);
		end

		if ischar(N)
			N=str2double(N);
		end		
	end



	dt=min([T/400 T/(10*N)]);
	if T/10000>dt
		if o~=1;
			dt=T/10000;
			disp(['Time step reduced to ' num2str(dt) '. See help to override.'])
			disp('Plot may look jagged when zoomed in.')
		else
			disp('Plot may look coarse. Time step size increased to prevent taking too long to run.')
		end
	end
	if N>500
		if o~=1;
			N=500;
			disp('Number of terms reduced to 500. See help to override.')
		else

			disp('This may take a long time to run unless you decrease N.')
		end
	end
	%dt=max([dt T/10000])
	t=0:dt:T*3;length(t);
	F=0*t+eval(a0)/2;
	for n=1:N
		if nargin>6
			if floor(n/2)==n/2
				a=aeven;
				b=beven;
			else
				a=aodd;
				b=bodd;
			end
			
		end
		F=F+eval(a)*cos(n*2*pi*t/T)+eval(b)*sin(n*2*pi*t/T);
	end

	plot(t,F)
	grid on

else
	vtb3_5(0,'-8/pi^2/n^2',0,0,0,20,10)
	title('Example 3.3.1')
	xlabel('time')
end

%Automatically check for updates
vtbchk

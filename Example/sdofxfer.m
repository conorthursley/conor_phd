	echo off
%	sdofxfer.m	plotting frequency responses of sdof model for different
%	damping values.  Calculates and plots magnitude and phase for a single
%	degree of freedom system over a range of damping values.
	
	clf;

	clear all;

%	assign values for mass, percentage of critical damping, and stiffnesses
%	zeta is a vector of damping values from 10% to 100% in steps of 10%

	m = 1;
	zeta = 0.1:0.1:1;			%   0.1 = 10% of critical
	k = 100;

	wn = sqrt(k/m);

%	Define a vector of frequencies to use, radians/sec.  The logspace command uses
%	the log10 value as limits, i.e. -1 is 10^-1 = 0.1 rad/sec, and 1 is
%	10^1 = 10 rad/sec.  The 400 defines 400 frequency points.

	w = logspace(-1,2,400);

%	pre-calculate the radians to degree conversion

	rad2deg = 180/pi;

%	define s as the imaginary operator times the radian frequency vector

	s = 1i*w;

%	define a for loop to cycle through all the damping values for calculating magnitude and phase

	for  cnt = 1:length(zeta)

%	define the frequency response to be evaluated

		xfer(cnt,:) = (1/m) ./ (s.^2 + 2*zeta(cnt)*wn*s + wn^2);

%	calculate the magnitude and phase of each frequency response

		mag(cnt,:) = abs(xfer(cnt,:));

		phs(cnt,:) = angle(xfer(cnt,:))*rad2deg;

	end

%	define a for loop to cycle through all the damping values for plotting magnitude

	for  cnt = 1:length(zeta)

	loglog(w,mag(cnt,:),'k-')
	title('SDOF frequency response magnitudes for zeta = 0.1 to 1.0 in steps of 0.1')
	xlabel('frequency, rad/sec')
	ylabel('magnitude')
	grid

	hold on

	end

	hold off
	
	grid on

	figure
%	define a for loop to cycle through all the damping values for plotting phase

	for  cnt = 1:length(zeta)

	semilogx(w,phs(cnt,:),'k-')
	title('SDOF frequency response phases for zeta = 0.1 to 1.0 in steps of 0.1')
	xlabel('frequency, rad/sec')
	ylabel('magnitude')
	grid

	hold on

	end

	hold off

	grid on
	
	


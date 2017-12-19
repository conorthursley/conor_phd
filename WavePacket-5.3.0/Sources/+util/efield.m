%--------------------------------------------------------------------------
%
% External electric field as a sequence of shaped pulses with
% constant or (linearly or quadratically) chirped carrier frequencies 
% with user-specified polarizations in the x/y plane and
% possibly with phase shifts.
%
%            N
% E(time) = sum f (time) * cos ( omega  * time + phase )
%           i=1  i                    i               i
%
% 
%                      /    (        ( time - delay_i )^2      )
%                     / exp ( -  ----------------------------- )
%                    /      (    2 ( fwhm_i / sqrt(8*ln2) )^2  )
% where f  (time) = <
%        i           \     2 (     time - delay_i    pi )
%                     \ sin  ( pi --------------- + --- )
%                      \     (        2 fwhm_i       2  )
%
% Note that the shifted sin^2 is actually rather a cos^2 but 
% we keep this notation for reasons of backward-compatibility
%
% Returning a structure with two vectors of field values (efield.x/y) 
% for any given vector (timesteps) of time values.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function efield = efield (timesteps)
global time

% Preallocate
efield.x = zeros(size(timesteps));
efield.y = zeros(size(timesteps));

% Loop over pulses
for ii = 1:time.efield.n_pulse

    % 0. Shift time wrt. pulse center
    tt0 = timesteps-time.efield.delay(ii);
    
    % 1. Calculate the envelope of the field
    switch lower(time.efield.shape(ii,:))
        case 'gauss'
            envelope = exp ( - tt0.^2 / ...
                ( 2 * (time.efield.fwhm(ii)/sqrt(8*log(2)))^2) );
        case 'recta'
            envelope = ones(size(tt0));
        case 'sin^2'
            envelope = cos ( pi/2 * tt0/time.efield.fwhm(ii) ).^2;
        case 'inter'
            envelope = interp1(time.efield.times{ii}, time.efield.fields{ii}, ...
                            tt0, time.efield.method, 0);
    end
    envelope = envelope * time.efield.ampli(ii);
    
    % 2. Calculate the oscillatory portion of the electric field.
    % constant or (linearly or quadratically) chirped frequencies
    frequency = time.efield.frequ(ii) * ones(size(tt0));
    frequency = frequency + time.efield.linear(ii) * tt0 ...
            + time.efield.quadratic(ii)/2 * tt0.^2;

    oscillat = cos ( frequency .* tt0 + time.efield.phase(ii) );
    
    % 3. Calculate the polarisation part of the electric field
    polar.x = cos(time.efield.polar(ii));
    polar.y = sin(time.efield.polar(ii));

    % 4. Put everything together
    where = find( abs(tt0) <= time.efield.length(ii)/2 );
    
    % Dressed state picture (Floquet) 
    if time.efield.dressed % Return half envelope
        efield.x(where) = efield.x(where) + ...
            0.5 * envelope(where) * polar.x;
        efield.y(where) = efield.y(where) + ...
            0.5 * envelope(where) * polar.y;
        
    % Bare state picture
    else % Return oscillating field
        efield.x(where) = efield.x(where) + ...
            envelope(where) .* oscillat(where) * polar.x;
        efield.y(where) = efield.y(where) + ...
            envelope(where) .* oscillat(where) * polar.y;
    end
    
end
end

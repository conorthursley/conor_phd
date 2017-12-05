%  Damped, driven harmonic oscillator: nonstiff ODE

function harmonic_sweep_example(fnumber)

    % Set up the initial condition
    omega = 1;      % natural frequency = sqrt(k/m)
    b = 0.3;        % drag coefficient s
    m = 1.0;        % mass?
    x0 = 1.0;       % initial position
    v0 = 1.0;       % initial velocity
    F0 = 1.0;       % strength of applied frequency?
    tBegin = 0;     % time begin
    tEnd = 80;      % time end

    % Applied frequency is chosen from random generator?
    rng('shuffle');
    omega_app = rand(1)*2.5;   % driving frequency


    % Use Runge-Kutta integrator to solve the ODE
    [t,w] = ode45(@derivatives, [tBegin tEnd], [x0 v0]);
    x = w(:,1);     % extract positions from first column of w matrix
    v = w(:,2);     % extract velocities from second column of w matrix

    xmax_norm = max(x)/x0   %normalize the displacement?
    omega_app_norm = omega_app/omega  %normalize the applied frequency


    % Write the outputs on a file?
    filenumber = num2str(fnumber);
    outfilename = sprintf ( '%s%s%s', 'XmaxFreq', filenumber, '.dat' );
    fileID = fopen(outfilename,'w');
    fprintf(fileID,'xmax_norm= %9.3f  omega_app_norm= %9.3f\n', xmax_norm, omega_app_norm);
    fclose(fileID);

    % Function to compute the  derivatives of dx/dt and dv/dt
    % The parameters m, b, F0, omega_app are from the main program?
    function derivs = derivatives(tf,wf)
        xf = wf(1);            % wf(1) stores x
        vf = wf(2);            % wf(2) stores v
        dxdt = vf;                                     % set dx/dt = velocity
        dvdt = - m*xf - b * vf + F0*cos(omega_app*tf);  % set dv/dt = acceleration
        derivs = [dxdt; dvdt];  % return the derivatives
    end

end
%% One dof example


%% Theoretical response for broadband excitation

mass1=0.1;		% [kg]
stiff1=1000;    % [N/m]
damp=0.0001;    % keep as a small number to fix solver errors


theory_freq=1:0.1:50;
theory_omega=2*pi*theory_freq;
theory_omega_0=sqrt(stiff1/mass1);

% say
eta=damp;
U_div_F = ones(size(theory_freq)) ./        ...
           ( (stiff1 - mass1*theory_omega.^2) +     ...
            1i*damp*theory_omega);

        
        
%% Load ansys harmonic results
[freq_ansys,ansys_harm_1dof_R,ansys_harm_1dof_I] = textread('ansys_1dof_ux.txt','%f%f%f');
ansys_harm_1dof_ux=ansys_harm_1dof_R+1i*ansys_harm_1dof_I;




%% Load ansys transient results
% This was from an analysis where there was white noise F=1N applied to the
% mass. The broadband excitation should give same results as the harmonic
% response
[time_trans,ansys_trans_ux] = textread('1dof_trans_ux2.txt','%f%f');
fs_trans=1/(time_trans(2)-time_trans(1));
FFTsize=1024;
[Pux_trans,Fux_trans]=pwelch(ansys_trans_ux,hanning(FFTsize),[],FFTsize,fs_trans);

[ansys_LP50Hz_force]=textread('filtered_white_data.txt','%f');
[Pfx_trans,Ffx_trans]=pwelch(ansys_LP50Hz_force,hanning(FFTsize),[],FFTsize,fs_trans);


[T_xfer,F_xfer] = tfestimate(ansys_LP50Hz_force,ansys_trans_ux,hanning(FFTsize),[],FFTsize,fs_trans);


%%% Load the ANSYS 10Hz force vector
load excitef.txt
excitef=excitef(:,4);
% calculate the PSD of the 10 Hz forcing function
[Pfx_trans10Hz,Ffx_trans10Hz]=pwelch(excitef,hanning(FFTsize),[],FFTsize,fs_trans);
% Calculate the xfer for the 10Hz forcing function
[T_xfer10Hz,F_xfer10Hz] = tfestimate(excitef,ansys_trans_ux,hanning(FFTsize),[],FFTsize,fs_trans);

figure
p_trans_10Hz=plot(Ffx_trans10Hz,10*log10(abs(Pfx_trans10Hz)),   ...
        F_xfer10Hz,10*log10(abs(T_xfer10Hz)),       ...
        freq_ansys,10*log10(abs(ansys_harm_1dof_ux)),'.',   ...
        Fux_trans,10*log10(abs(Pux_trans))  );



%% Theoretical response for 10Hz excitation
% see
% https://en.wikipedia.org/wiki/Harmonic_oscillator
% Sinusoidal driving force
omega_tone=2*pi*10;
Z_m=sqrt(  (2*theory_omega_0*eta)^2 + 1/omega_tone^2*(theory_omega_0^2-omega_tone^2)^2);
theta=atan( (2*omega_tone*theory_omega_0*eta)/(omega_tone^2-theory_omega_0^2) );
time=time_trans;
x_tone_theory = 1 / (mass1*Z_m*omega_tone) * sin(omega_tone*time+theta);

%%% Calculate the PSD of the response
[Ptheory_tone,Ftheory_tone]=pwelch(x_tone_theory,hanning(FFTsize),[],FFTsize,fs_trans);


%% Plot results

figure
p_theory=plot(theory_freq,10*log10(abs(U_div_F)));
xlabel('Frequency (Hz)');
ylabel('dB)')
%legend('Theory')

hold on
p_theory_tone=plot(Ftheory_tone,10*log10(abs(Ptheory_tone)));

hold on
p_ansys=plot(freq_ansys,10*log10(abs(ansys_harm_1dof_ux)),'.');

hold on
p_trans=plot(Fux_trans,10*log10(abs(Pux_trans)),'--',   ...
            Ffx_trans,10*log10(abs(Pfx_trans)), ...
            F_xfer,10*log10(abs(T_xfer))    );
legend( 'Theory Harmonic',   ...
        'Ansys Harmonic',     ...
        'Ansys Trans UX',   ...
        'Ansys Trans FX',   ...
        'Ansys Trans xfer');

axis([0 50 -120 0])

%% Script to plot transfer function from the chain.
% Carl Howard 26/1/2018


%% Load the results from transient analysis
[time1,disp_n1] = textread('disp_node2.txt','%f%f');
[time2,disp_n_end] = textread('disp_node_end.txt','%f%f');
[time3,ignore1,ignore2,force] = textread('force.txt','%f%f%f%f');
clear ignore1 ignore2


excitation=disp_n1;
response=disp_n_end;

fs=1./(time1(3)-time1(2));

%[Txy,F3] = tfestimate(force,signal(:,2),hanning(length(signal)),[],FFTsize,fs)
FFTsize=1024;
[Txy,Fxy] = tfestimate(excitation,response,hanning(FFTsize),[],FFTsize,fs);
%[Txy,F3] = tfestimate(excitation,response,512,[],[],fs);


[Pxx,Fxx]=pwelch(excitation,hanning(FFTsize),[],FFTsize,fs);
[Pyy,Fyy]=pwelch(response,hanning(FFTsize),[],FFTsize,fs);
[Pzz,Fzz]=pwelch(force,hanning(FFTsize),[],FFTsize,fs);

%% Load results from harmonic analysis
[freq1,harm_disp_ux1_R,harm_disp_ux1_I] = textread('ansys_harmonic_ux1.txt','%f%f%f');
[freq2,harm_disp_uxn_R,harm_disp_uxn_I] = textread('ansys_harmonic_uxn.txt','%f%f%f');
[freq3,harm_force_R,harm_force_I] = textread('ansys_harmonic_force.txt','%f%f%f');
[freq4,harm_xfer_R,harm_xfer_I] = textread('ansys_harmonic_xfer.txt','%f%f%f');

%% Convert into complex variables and delete originals
harm_disp_ux1=harm_disp_ux1_R+i*harm_disp_ux1_I;
clear harm_disp_ux1_R harm_disp_ux1_I
harm_disp_uxn=harm_disp_uxn_R+i*harm_disp_uxn_I;
clear harm_disp_uxn_R harm_disp_uxn_I
harm_force=harm_force_R+i*harm_force_I;
clear harm_force_R harm_force_I
harm_xfer=harm_xfer_R+i*harm_xfer_I;
clear harm_xfer_R harm_xfer_I



%% This is used to test the FFT sine wave is working properly
figure
testsine=1*sin(2*pi*50*time1);
[Psine,Fsine]=pwelch(testsine,hanning(FFTsize),[],FFTsize,fs);
plot(Fsine,10*log10(Psine));
xlabel('Frequency');
ylabel('(dB)')
legend('pwlech 50Hz sine wave')


%% Plot the transient results
figure
p_trans=plot(Fxx,10*log10(Pxx),'-',Fyy,10*log10(Pyy),'-',Fzz,10*log10(Pzz),'-',Fxy,10*log10(abs(Txy)),'-');
xlabel('Frequency');
ylabel('(dB)')
legend('PSD disp excitation','PSD disp response','force','disp end node / disp 1st node')

%% Plot the ANSYS harmonic results
hold on
p10=plot(   freq1,20*log10(abs(harm_disp_ux1)),'-.',    ...
            freq2,20*log10(abs(harm_disp_uxn)),'-.', ...
            freq3,20*log10(abs(harm_xfer)),'-.', ...
            freq4,20*log10(abs(harm_force)),'-.'    );

axis([0 100 -200 20])

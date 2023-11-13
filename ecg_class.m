%The following MATLAB code is used to plot
%the ECG signal and its Periodogram spectrum.
clc
clear
%Load Sample ECG Data downloaded from the web site
%https://www.physionet.org/physiobank/database/mitdb/
load MITBIH_ECG.mat
%%
Fs = 320; % Sampling frequency
T = 1/Fs; % Sampling period
L = length(ECGN(:,1)); % Length of signal
t = (0:L-1)*T; % Time vector
%%
%Plot NORMAL ECG Signal
subplot(2,1,1);
plot(ECGN(:,1));
title(’Original ECG Signal’)
xlabel(’Number of Samples’)
ylabel(’Amplitude’)
ylabel(’X(t)’)
axis ([0 L –inf inf])
%%
% Obtain the periodogram of the signal using periodogram. Plot the periodogram.
[Pxx,F] = periodogram(ECGN(:,1),hamming(L),L,Fs);
subplot(2,1,2);
f = Fs*(0:(L/2))/L;
plot(f,Pxx)
title(’Periodogram of ECG Signal’)
xlabel(’Frequency (Hz)’);
ylabel(’jP(f)j’)
axis ([0 50 –inf inf])


noverlap = L/4;
order = 34;
%%
% Obtain the Yule-Walker AR Spectrum of the Normal ECG signal using pyulear.
[Pxx,F] = pyulear(ECGN(:,1),order,segmentLength,Fs);
% Plot the Yule-Walker AR Spectrum of the signal
subplot(2,2,1)
plot(F,10*log10(Pxx))
ylabel(’PSD (dB/Hz)’);
xlabel(’Frequency (Hz)’);
legend(’Normal ECG’)
axis([0 L/2 –inf inf])
%title(’Yule-Walker AR Spectrum of ECG Signals’)
%%
% Obtain the Yule-Walker AR Spectrum of the ECG signal with APC using pyulear.
[Pxx,F] = pyulear(ECGAPC(:,1),order,segmentLength,Fs);
% Plot the Yule-Walker AR Spectrum of the signal
subplot(2,2,2)
plot(F,10*log10(Pxx))
ylabel(’PSD (dB/Hz)’);
xlabel(’Frequency (Hz)’);
legend(’APC’)
axis([0 L/2 –inf inf])
%title(’Yule-Walker AR’)
%%
% Obtain the Yule-Walker AR Spectrum of the ECG signal with LBBB using pyulear.
[Pxx,F] = pyulear(ECGLBBB(:,1),order,segmentLength,Fs);
% Plot the Yule-Walker AR Spectrum of the signal
subplot(2,2,3)
plot(F,10*log10(Pxx))
ylabel(’PSD (dB/Hz)’);
xlabel(’Frequency (Hz)’);
legend(’LBBB’)
axis([0 L/2 –inf inf])
%title(’Yule-Walker AR’)
%%
% Obtain the Yule-Walker AR Spectrum of the ECG signal with PVC using pyulear.
[Pxx,F] = pyulear(ECGPVC(:,1),order,segmentLength,Fs);
% Plot the Yule-Walker AR Spectrum of the signal
subplot(2,2,4)
plot(F,10*log10(Pxx))
ylabel(’PSD (dB/Hz)’);
xlabel(’Frequency (Hz)’);
legend(’PVC’)
axis([0 L/2 –inf inf])
%title(’Yule-Walker AR’)

%%
% Wigner-Ville time-frequency distribution
figure
[tfr,t,f] = tfrwv(ECGN(:,1));
mesh(t, f,abs(tfr(:,:)))
ylabel(’Frequency(Hz) ’); xlabel(’Time(msec)’);zlabel(’dB’);
title(’Wigner Ville of ECG Signal’)
%%
%Pseudo Wigner-Ville time-frequency distribution.
figure
[tfr,t,f] = tfrpwv(ECGN(:,1));
mesh(t, f,abs(tfr(:,:)))
ylabel(’Frequency(Hz) ’); xlabel(’Time(sec)’);zlabel(’dB’);
%axis([-inf inf 0 100 –inf inf])
title(’Pseudo Wigner Ville of ECG Signal’)
%%
%Smoothed Pseudo Wigner-Ville time-frequency distribution.
figure
[tfr,t,f] = tfrspwv(ECGN(:,1));
mesh(abs(tfr))
ylabel(’Frequency(Hz) ’); xlabel(’Time(sec)’);zlabel(’dB’);
title(’Smoothed Pseudo Wigner Ville of ECG Signal’)

% Wigner Ville
figure
[tfr,t,f] = tfrcw(ECGN(:,1));
mesh(t, f,abs(tfr(:,:)))
ylabel(’Frequency(Hz) ’); xlabel(’Time(msec)’);zlabel(’dB’);
title(’Choi-Williams time-frequency distribution of ECG Signal’)


% Obtain and plot the CWT of Normal ECG Signal using the analytic Morlet wavelet.
t = 0:DT:(numel(ECGN(:,1))*DT)–DT;
f0 = 5/(2*pi);
scales = helperCWTTimeFreqVector(20,500,f0,0.001,32);
cwtecg = cwtft({ECGN(:,1),0.001},’wavelet’,’morl’,’scales’,scales);
subplot(2,2,1);
helperCWTTimeFreqPlot(cwtecg.cfs,ECGN(:,1),cwtecg.frequencies,...
’surf’,’CWT of Normal ECG Signal’,’Time (secs)’,’Frequency (Hz)’)
colorbar off
%%
% Obtain and plot the CWT of APC ECG Signal using the analytic Morlet wavelet.
t = 0:DT:(numel(ECGAPC(:,1))*DT)–DT;
f0 = 5/(2*pi);
scales = helperCWTTimeFreqVector(20,500,f0,0.001,32);
cwtecg = cwtft({ECGAPC(:,1),0.001},’wavelet’,’morl’,’scales’,scales);
subplot(2,2,2);
helperCWTTimeFreqPlot(cwtecg.cfs,ECGAPC(:,1),cwtecg.frequencies,...
’surf’,’CWT of APC ECG Signal’,’Time (secs)’,’Frequency (Hz)’)
colorbar off
%%
% Obtain and plot the CWT of PVC ECG Signal using the analytic Morlet wavelet.
t = 0:DT:(numel(ECGPVC(:,1))*DT)–DT;
f0 = 5/(2*pi);
scales = helperCWTTimeFreqVector(20,500,f0,0.001,32);
cwtecg = cwtft({ECGPVC(:,1),0.001},’wavelet’,’morl’,’scales’,scales);
subplot(2,2,3);
helperCWTTimeFreqPlot(cwtecg.cfs,ECGPVC(:,1),cwtecg.frequencies,...
’surf’,’CWT of PVC ECG Signal’,’Time (secs)’,’Frequency (Hz)’)
colorbar off
%%
% Obtain and plot the CWT of RBBB ECG Signal using the analytic Morlet wavelet.
t = 0:DT:(numel(ECGRBBB(:,1))*DT)–DT;
f0 = 5/(2*pi);
scales = helperCWTTimeFreqVector(20,500,f0,0.001,32);
cwtecg = cwtft({ECGRBBB(:,1),0.001},’wavelet’,’morl’,’scales’,scales);
subplot(2,2,4);
helperCWTTimeFreqPlot(cwtecg.cfs,ECGRBBB(:,1),cwtecg.frequencies,...
’surf’,’CWT of RBBB ECG Signal’,’Time (secs)’,’Frequency (Hz)’)
colorbar off
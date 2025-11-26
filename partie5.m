clear; clc; close all;
Fs = 8000; 
load('fcno03fz.mat');

% Initialisation des variables
s = fcno03fz'; %signal non bruité
N = length(s);
t = (0:N-1) / Fs;
RSB = 10; 

% Bruitage du signal
P_signal=mean(s.^2);
P_bruit_ = P_signal / (10^(RSB / 10));
P_bruit=sqrt(P_bruit_); %ecart-type
bruit=P_bruit*randn(size(s));
signal_bruite=s+bruit;

% Découpage en trames
n_s = round(0.03 * Fs);   % trames de 30 ms
pas = round(n_s / 2);     % recouvrement 50 %
NFFT = 2*n_s-1;
N_welch = floor(n_s/10);

nbr_trames = floor((length(signal_bruite)-n_s)/pas)+1;
if length(signal_bruite) < (nbr_trames-1)*pas + n_s
    signal_bruite = [signal_bruite zeros(1, (nbr_trames-1)*pas + n_s - length(s))];
end

w = hann(n_s)';

x_trames_origin = zeros(nbr_trames, n_s);
x_trames_bruite = zeros(nbr_trames, n_s);
x_trames_rehausse = zeros(nbr_trames, n_s);

spectres_puissances_origin = zeros (nbr_trames, NFFT);
spectres_puissances_bruite = zeros (nbr_trames, NFFT);
spectres_puissances_rehausse = zeros (nbr_trames, NFFT);

%Estimation de la DSP du bruit
DSP_welch = Perio_Welch(bruit, NFFT);

for i=1:nbr_trames
    % décomposition et fenetrage
    x_trames_bruite(i, :) = signal_bruite((i-1)*pas+1 : (i-1)*pas + n_s) .* w;
    x_trames_origin(i, :) = s((i-1)*pas+1 : (i-1)*pas + n_s) .* w;

    %calcul des spectres de puissance (original et bruité)
    spectres_puissances_origin(i,:) = fft(x_trames_origin(i,:), NFFT);
    spectres_puissances_origin(i,:) = (1/NFFT)*abs(spectres_puissances_origin(i,:)).^2;

    spectres_puissances_bruite(i,:) = fft(x_trames_bruite(i,:), NFFT);
    spectres_puissances_bruite(i,:) = (1/NFFT)*abs(spectres_puissances_bruite(i,:)).^2;

    %Spectre réhaussé
    for j = 1:NFFT
        spectres_puissances_rehausse(i,j) = spectres_puissances_bruite(i,j) - DSP_welch(j);
        if (spectres_puissances_rehausse(i,j) < 0)
            spectres_puissances_rehausse(i,j) = 0;
        end
    end
    phase_bruite = angle(fft(x_trames_bruite(i,:), NFFT)); %phase du signal bruite
    TF_signal_rehausse = sqrt(NFFT * spectres_puissances_rehausse(i,:)) .* exp(1j * phase_bruite); %TF du signal rehausse

    %Calcul du signal temporel rehausse
    x_temp = real(ifft(TF_signal_rehausse, NFFT));
    x_trames_rehausse(i,:) = x_temp(1:n_s);
end

% Représentations
time_axis = (1/Fs)*(1:n_s);
freq_axis = -Fs/2:Fs/NFFT:Fs/2-Fs/NFFT;

figure;
plot(time_axis, x_trames_origin(floor(nbr_trames/2),:), "red");
hold on;
plot(time_axis, x_trames_bruite(floor(nbr_trames/2),:), "green");
plot(time_axis, x_trames_rehausse(floor(nbr_trames/2),:), "blue");
hold off;
title('Représentation des signaux temporels');
legend('Signal original', 'Signal bruité', 'Signal rehaussé');


figure;
plot(freq_axis, fftshift(spectres_puissances_origin(floor(nbr_trames/2),:)), "red");
hold on;
plot(freq_axis, fftshift(spectres_puissances_bruite(floor(nbr_trames/2),:)), "green");
plot(freq_axis, fftshift(spectres_puissances_rehausse(floor(nbr_trames/2),:)), "blue");
hold off;
title('Représentation des spectres de puissance');
legend('Signal original', 'Signal bruité', 'Signal rehaussé');


%addition recouvrement
y = zeros(1,length(s));
w_sum = zeros(1,length(s));

for i = 1:nbr_trames
    idx = (i-1)*pas+1 : (i-1)*pas + n_s;
    y(idx) = y(idx) + x_trames_rehausse(i,:);
    w_sum(idx) = w_sum(idx) + w;       % fenêtre de Hann
end

eps_val = 1e-12;                
w_sum(w_sum < eps_val) = eps_val;
for j = 1:length(s)
    if (w_sum(j) ~= eps_val)
        y(j) = y(j) / w_sum(j);
    end
end
% y = y ./ w_sum;
for i = 1:nbr_trames
    if (abs(y(i)) > 30000)
        y(i) = 0;
    end
end
y(~isfinite(y)) = 0;

% Paramètres pour le Spectrogramme
win_len = round(0.030 * Fs); % 30ms
overlap = round(0.5 * win_len); % 50% recouvrement
nfft = 1024; % Nombre de points FFT

%Représentations
figure;

subplot(3,1,1);
plot(t, y);
xlabel('temps (s)');
ylabel('Amplitude');
title("Signal temporel réhaussé");

subplot(3,1,2);
spectrogram(y, hamming(win_len), overlap, nfft, Fs, 'yaxis');
title('Spectrogramme');

subplot(3,1,3);
plot(t, s-y);
xlabel('temps (s)');
ylabel('Amplitude');
title("diff");







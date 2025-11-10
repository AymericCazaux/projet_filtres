
Fs = 8000; 
load('fcno03fz.mat');

s = fcno03fz;
%sound(s, Fs);
N = length(s);
t = (0:N-1) / Fs;

RSB_test = 10; 

s_bruite = bruiter(s, RSB_test);


% --- Paramètres pour le Spectrogramme ---
win_len = round(0.030 * Fs); % 30ms
overlap = round(0.5 * win_len); % 50% recouvrement
nfft = 1024; % Nombre de points FFT

% --- Figure 1: Signal Original ---
figure(1);
sgtitle('Signal de Parole Original');

% Plot temporel
subplot(2, 1, 1);
plot(t, s);
title('Représentation Temporelle');
xlabel('Temps (s)');
ylabel('Amplitude');
xlim([t(1) t(end)]); 

% Spectrogramme
subplot(2, 1, 2);
spectrogram(s, hamming(win_len), overlap, nfft, Fs, 'yaxis');
title('Spectrogramme');

% --- Figure 2: Signal Bruité ---
figure(2);
sgtitle(['Signal de Parole Bruité (RSB = 10 dB)']);

% Plot temporel
subplot(2, 1, 1);
plot(t, s_bruite);
title('Représentation Temporelle');
xlabel('Temps (s)');
ylabel('Amplitude');
xlim([t(1) t(end)]); 

% Spectrogramme
subplot(2, 1, 2);
spectrogram(s_bruite, hamming(win_len), overlap, nfft, Fs, 'yaxis');
title('Spectrogramme');


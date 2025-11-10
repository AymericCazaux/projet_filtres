%La réponse h(k) a une durée finie car on n'a que deux impulsions à deux moments distincts. C'est donc un filtre à Réponse Impulsionnelle Finie (RIF)
%La sortie ne dépend que des entrées présentes et passées donc il est causal.
%IL est stable car un filtre RIF est toujours stable

%Fonction de Transfert (H(z)) : C'est la Transformée en Z de h(k).Donc : H(z)=1+z−k0

%

Fs=8000;
load('fcno03fz.mat');

s = fcno03fz;
k0=10;
p=[1];
num=[1 zeros(1, k0-1), 1];
filtre=filter(num,p, s);

% --- Paramètres pour le Spectrogramme ---
win_len = round(0.030 * Fs); % 30ms
overlap = round(0.5 * win_len); % 50% recouvrement
nfft = 1024; % Nombre de points FFT

% --- Figure 1: Signal Original ---
N = length(s);
t = (0:N-1) / Fs;
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

% --- Figure 2: Signal Filtré ---
figure(2);
sgtitle(['Signal de parole filtré']);

% Plot temporel 
subplot(2, 1, 1);
plot(t, filtre);
title('Représentation Temporelle');
xlabel('Temps (s)');
ylabel('Amplitude');
xlim([t(1) t(end)]); 

% Spectrogramme
subplot(2, 1, 2);
spectrogram(filtre, hamming(win_len), overlap, nfft, Fs, 'yaxis');
title('Spectrogramme');

%il faut bien enlever les zones mortes dans les figures et aligner temporel
%et spectogramme (enlever si besoin la barre)(utiliser property inspection
%peut-être.

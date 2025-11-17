 clc; clear; close all;

%Variables
N = 1000; %nombre de points
Tech = 10^(-3);  %période d'échantillonage
fech = 1/Tech;
m = 0;  %moyenne du bruit blanc gaussien
sigma2 = 1;  %variance
time_axis = -(N-1)*Tech:Tech:(N-1)*Tech;
freq_axis = -fech/2+fech/(2*N):fech/(2*N):fech/2-fech/(2*N);


%Génération du bruit blanc gaussien centré
bruit = m + sqrt(sigma2)*randn(1,N);

%% fonction d'auto-correlation théorique
Re_bruit_th = zeros(1, 2*N -1);
Re_bruit_th(N+1) = sigma2;

%fonction d'auto-correlation estimée

Re_non_biaisee = xcorr(bruit, 'unbiased');
Re_biaisee = xcorr(bruit, 'biased');


%représentation des fonction d'autocorrelation
figure;
subplot(3,1,1);
plot(time_axis, Re_bruit_th);
title('Fonction théorique');
grid on;

subplot(3,1,2);freq_axis = -fech/2+fech/(2*N):fech/(2*N):fech/2-fech/(2*N);
plot(time_axis, Re_biaisee);
title('Estimateur biaisé');
grid on;

subplot(3,1,3);
plot(time_axis, Re_non_biaisee);
title('Estimateur non biaisé');
grid on;

%% DSP théorique
DSP_th = sigma2* ones(1,2*N -1);

%spectres de puissance
spectre_puissance = fft(bruit, 2*N-1);
spectre_puissance = (1/(2*N-1))*abs(spectre_puissance).^2;
spectre_puissance = fftshift(spectre_puissance);

%représentation des spectres et dsp
figure;
subplot(2,1,1);
plot(freq_axis, DSP_th);
title('DSP théorique');

subplot(2,1,2);
plot(freq_axis, spectre_puissance);
title('spectre de puissance');

%% calcul des periodogrammes et du correlogramme
NFFT = 100;
Daniell = Perio_Daniell(bruit, N, 20);
Barlett = Perio_Barlett(bruit, NFFT);
Welch = Perio_Welch(bruit, NFFT);
Corello = Correlogramme(Re_biaisee, N);

freq_axis_correlo = -fech/2:fech/N:fech/2-fech/N;
freq_axis_perio = -fech/2:fech/NFFT:fech/2-fech/NFFT;


%Représentation de ces derniers
figure;
subplot(4,1,1);
plot(freq_axis_correlo, fftshift(Daniell));
title('Periodogramme de Daniell')

subplot(4,1,2);
plot(freq_axis_perio, fftshift(Barlett));
title('Periodogramme de Barlett');

subplot(4,1,3);
plot(freq_axis_perio, fftshift(Welch));
title('Periodogramme de Welch');

subplot(4,1,4);
plot(freq_axis_correlo, fftshift(Corello));
title('Correlogramme');


%% Platitude spectrale

nbr_reps = 10;
platitude_spectrale = zeros(1, nbr_reps);
platitude_spectrale_welch = zeros(1, nbr_reps);
platitude_spectrale_theorique = zeros(1, nbr_reps);
platitude_spectrale_ordrep = zeros(1, nbr_reps);

for j = 1:nbr_reps

    %Génération du bruit blanc gaussien centré
    bruit = m + sqrt(sigma2)*randn(1,N);
    %spectres de puissance
    spectre_puissance = fft(bruit, 2*N-1);
    spectre_puissance = (1/(2*N-1))*abs(spectre_puissance).^2;
    spectre_puissance = fftshift(spectre_puissance);
    %Periodogramme de Welch
    Welch = Perio_Welch(bruit, NFFT);

    %Moyenne arithmétique (spectre)
    moy_arith_s = 0;
    for i = 1:length(spectre_puissance)
        moy_arith_s = moy_arith_s + spectre_puissance(i);
    end
    moy_arith_s = moy_arith_s / length(spectre_puissance);

    %Moyenne arithmétique (welch)
    moy_arith_w = 0;
    for i = 1:length(Welch)
        moy_arith_w = moy_arith_w + Welch(i);
    end
    moy_arith_w = moy_arith_w / length(Welch);
    
    %Moyenne géométrique (spectre)
    moy_geo_s = 1;
    for i = 1:length(spectre_puissance)
        moy_geo_s = moy_geo_s * spectre_puissance(i);
    end
    moy_geo_s = moy_geo_s^(1/length(spectre_puissance));

    %Moyenne géométrique (welch)
    moy_geo_w = 1;
    for i = 1:length(Welch)
        moy_geo_w = moy_geo_w * Welch(i);
    end
    moy_geo_w = moy_geo_w^(1/length(Welch));

    %Moyenne d ordre p
    p = 2;
    moy_p = 0;
    for i =1:length(spectre_puissance)
        moy_p = moy_p + spectre_puissance(i)^p;
    end
    moy_p = moy_p / length(spectre_puissance);
    moy_p = moy_p^(1/p);
    
    platitude_spectrale(j) = moy_geo_s / moy_arith_s;
    platitude_spectrale_ordrep(j) = moy_p / moy_arith_s;
    platitude_spectrale_welch(j) = moy_geo_w / moy_arith_w;
end


moy_plat_spectrale = 0;
moy_plat_spectrale_welch = 0;
moy_plat_spectrale_p = 0;

for j = 1:nbr_reps
    moy_plat_spectrale = moy_plat_spectrale + platitude_spectrale(j)/nbr_reps;
    moy_plat_spectrale_welch = moy_plat_spectrale_welch + platitude_spectrale_welch(j)/nbr_reps;
    moy_plat_spectrale_p = moy_plat_spectrale_p + platitude_spectrale_ordrep(j)/nbr_reps;
end

moy_plat_spectrale_p
moy_plat_spectrale
moy_plat_spectrale_welch


    

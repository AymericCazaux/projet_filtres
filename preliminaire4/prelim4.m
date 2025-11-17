clear; clc; close all;

%%Procédure d'addition recouvrement
load('fcno03fz.mat');

%Initialisation des variables

s = fcno03fz(10000:20000)';
pas = floor(length(s)/10);
n_s = floor(2*pas); %recouvrement à 50%
len_bord = n_s-pas;

y = addition_recouvrement(s, pas, n_s);

z = s-y;

% z = abs(z);
% err_moy = mean(z)

figure;
subplot(3,1,1);
plot(s);

subplot(3,1,2);
plot(y);

subplot(3,1,3);
plot(z);


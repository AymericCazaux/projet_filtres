function [signal_bruite]=bruiter(signal, RSB)
    P_signal=mean(signal.^2);
    P_bruit_ = P_signal / (10^(RSB / 10));
    P_bruit=sqrt(P_bruit_); %ecart-type
    bruit=P_bruit*randn(size(signal));
    signal_bruite=signal+bruit;
end


%Spectrogramme Original : L'allure est cohérente d'un signal de parole. On peut clairement faire la différence entre périodes de silence (rès basse énergie) et périodes de parole active. On y distingue les formants (les bandes horizontales jaunes à haute énergie ) ainsi que la structure harmonique (les fines stries verticales) dans les sons voisé
%Spectrogramme Bruité (10 dB) : L'ajout du bruit blanc, qui est stationnaire et a une densité spectrale de puissance théoriquement plate, se traduit par un "tapis" de bruit (le fond vert/jaune) à un niveau d'énergie quasi constant sur toutes les fréquences et sur toute la durée du signal.

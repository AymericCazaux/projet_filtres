function y = Perio_Daniell(x, NFFT, m)
    
    %Initialisation  
    fft_x = fft(x, NFFT);
    Perio_x = (abs(fft_x).^2) /NFFT;
    y = Perio_x;

    %Lissage
    lissage = 0;
    compteur = 0;
    for i =1:length(y)
        for j= -m:m
            if ((i+j) >= 1 && (i+j) <= length(y))
                lissage = lissage + y(i+j);
                compteur = compteur + 1;
            end
        end
        lissage = lissage / compteur;
        y(i) = lissage;
        lissage = 0;
        compteur = 0;
    end
end
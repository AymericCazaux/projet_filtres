function y = addition_recouvrement(s, pas, n_s)

    % s: signal d'origine
    % pas: pas entre chaque trame
    % n_s: longueur d'une trame
    
    s = s(:)';
    nbr_trames = floor((length(s)-n_s)/pas)+1;
    if length(s) < (nbr_trames-1)*pas + n_s
        s = [s zeros(1, (nbr_trames-1)*pas + n_s - length(s))];
    end
    w = hann(n_s)';
    x_trames = zeros(nbr_trames, n_s);

    for i=1:nbr_trames
        % décomposition et fenetrage
        x_trames(i, :) = s((i-1)*pas+1 : (i-1)*pas + n_s) .* w;
    end

    % Reconstruction du signal
    y = zeros (1, length(s));
    for i=1:nbr_trames
        y((i-1)*pas+1 : (i-1)*pas + n_s) = y((i-1)*pas+1 : (i-1)*pas + n_s) + x_trames(i, :);
    end
end




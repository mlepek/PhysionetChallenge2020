function [newann, filling] = holefilling(oldann, fillalg, Fs, tolerancja,signal)
% Przeznaczenie:
% funkcja wype³niaj¹ca przerwy w anotacjach w metodach hierarchicznych.
% argumenty wejœciowe:
% oldann - istniej¹ce wczeœniej anotacje
% fillalg - algorytm wype³niajacy przerwy
% Fs - czeœtotliwoœæ próbkowania
% tolerancja - czas refrakcji
% signal - sygna³
% Funkcja zwraca:
% newann - indeksy anotacji
% uwagi:
% sygna³y i anotacje musz¹ byæ zgodne z formatem u¿ywanym w WFDB Toolbox
% Autor: Antonina Pater, 2018
% Ostatnia modyfikacja: 1.01.19

oldann = [1; oldann; length(signal)];
RR = diff(oldann);
% znaleœæ ró¿nicê o min 1400 ms,
holeinx = find(RR >= 1.4*Fs);
% policzyæ dla tego wycinka 2 algorytm
filling = [];
for ii = 1:length(holeinx)
    fillingtmp = fillalg(fillalg>oldann(holeinx(ii))...
        &fillalg<RR(holeinx(ii))+oldann(holeinx(ii)));
    % usun¹æ wype³nienia za blisko wa¿niejszego detektora
    for a = 1:length(fillingtmp)
        if (fillingtmp(a) < (RR(holeinx(ii)) + Fs*tolerancja)) ...
                | (fillingtmp(a) > (RR(holeinx(ii)) + ... 
                oldann(holeinx(ii)) - Fs*tolerancja))
            fillingtmp(a) = 0;
        end
    end

    fillingtmp(fillingtmp==0) = [];
    filling = [filling; fillingtmp];
end

oldann = [oldann; filling];
newann = sort(oldann);
newann(newann==1) = [];
newann(newann >= length(signal)) = [];
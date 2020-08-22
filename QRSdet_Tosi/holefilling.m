function [newann, filling] = holefilling(oldann, fillalg, Fs, tolerancja,signal)
% Przeznaczenie:
% funkcja wype�niaj�ca przerwy w anotacjach w metodach hierarchicznych.
% argumenty wej�ciowe:
% oldann - istniej�ce wcze�niej anotacje
% fillalg - algorytm wype�niajacy przerwy
% Fs - cze�totliwo�� pr�bkowania
% tolerancja - czas refrakcji
% signal - sygna�
% Funkcja zwraca:
% newann - indeksy anotacji
% uwagi:
% sygna�y i anotacje musz� by� zgodne z formatem u�ywanym w WFDB Toolbox
% Autor: Antonina Pater, 2018
% Ostatnia modyfikacja: 1.01.19

oldann = [1; oldann; length(signal)];
RR = diff(oldann);
% znale�� r�nic� o min 1400 ms,
holeinx = find(RR >= 1.4*Fs);
% policzy�dla tego wycinka 2 algorytm
filling = [];
for ii = 1:length(holeinx)
    fillingtmp = fillalg(fillalg>oldann(holeinx(ii))...
        &fillalg<RR(holeinx(ii))+oldann(holeinx(ii)));
    % usun�� wype�nienia za blisko wa�niejszego detektora
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
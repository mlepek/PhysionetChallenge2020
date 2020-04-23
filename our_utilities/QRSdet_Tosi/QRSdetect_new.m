function [newann, newL]  = QRSdetect_new(RR, h, L, Fs, ann, tolerancja, signal)
% Przeznaczenie:
% funkcja jest wywo�ywana w waveletdet.m, s�u�y do uzupe�niania anotacji
% przez kolejny algorytm
% argumenty wej�ciowe:
% RR - d�ugi�ci interwa��w
% h - warto�c iloczynu h
% L - warto�� progu L
% Fs - cze�totliwo�� pr�kowania sygna�u
% ann - anotacje
% tolerancja - zakres w kt�rym anotacje s� uznawane za wskazuj�ce na ten
% sam zesp� QRS
% signal - sygna� EKG
% uwagi:
% do poprawnego dzia�ania funkcji konieczne jest zainstalowanie biblioteki
% WFDB Toolbox
% sygna�y i anotacje musz� by� zgodne z formatem u�ywanym w WFDB Toolbox
% Autor: Antonina Pater, 2018
% Ostatnia modyfikacja: 1.01.19

newL = 0.5 * L; % pr�g mniejszy o po�ow�
suspRRinx = find(RR >= 2*Fs); %podejrzanie d�ugie interwa�y
newann = ann;

for j = 1:length(suspRRinx)
    
    anntmp = find(h(ann(suspRRinx(j)):ann(suspRRinx(j))+RR(suspRRinx(j))) >= newL);
    newann = [newann; ann(suspRRinx(j)) + anntmp];
end

%przejscie od prawej do lewej
for j = length(suspRRinx):-1:1
    if suspRRinx(j)>1
    anntmp = find(h(ann(suspRRinx(j))+RR(suspRRinx(j)):ann(suspRRinx(j)))>= newL);
    newann = [newann; ann(suspRRinx(j)') + anntmp];
    end
end

sort(newann);
newann = unique(newann);

newann(newann>length(signal)) = [];
% petla po wszystkich anotacjach i zostawienie tylko tych, ktore wskazuja
% na lokalne maksimum.
    signalAtAnn = signal(newann);
    newannToremove = [];
 for j = 1:length(newann)
     newannTemp = newann(abs(newann(j)-newann)<tolerancja*Fs);
     signalAtAnnTemp = signal(newann(j));
     if(signalAtAnnTemp <  max(signal(newannTemp)))
         newannToremove = [newannToremove ; j];
     end
 

 end
 
 newann(newannToremove) = []; 
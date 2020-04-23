function [newann, newL]  = QRSdetect_new(RR, h, L, Fs, ann, tolerancja, signal)
% Przeznaczenie:
% funkcja jest wywo³ywana w waveletdet.m, s³u¿y do uzupe³niania anotacji
% przez kolejny algorytm
% argumenty wejœciowe:
% RR - d³ugiœci interwa³ów
% h - wartoœc iloczynu h
% L - wartoœæ progu L
% Fs - czeœtotliwoœæ prókowania sygna³u
% ann - anotacje
% tolerancja - zakres w którym anotacje s¹ uznawane za wskazuj¹ce na ten
% sam zespó³ QRS
% signal - sygna³ EKG
% uwagi:
% do poprawnego dzia³ania funkcji konieczne jest zainstalowanie biblioteki
% WFDB Toolbox
% sygna³y i anotacje musz¹ byæ zgodne z formatem u¿ywanym w WFDB Toolbox
% Autor: Antonina Pater, 2018
% Ostatnia modyfikacja: 1.01.19

newL = 0.5 * L; % próg mniejszy o po³owê
suspRRinx = find(RR >= 2*Fs); %podejrzanie d³ugie interwa³y
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
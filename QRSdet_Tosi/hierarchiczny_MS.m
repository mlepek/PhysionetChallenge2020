function newann = hierarchiczny_MS(signal, Fs,ECG_RawData,HRVparams)
% Przeznaczenie:
% hierarchiczna metoda decyzyjna ustalaj¹ca hierarchiê algorytmów
% podstawowych wg ustalonej hierarchii
% argumenty wejœciowe:
% name - nazwa sygna³u
% Funkcja zwraca:
% newann - indeksy anotacji
% U¿ywane funkcje:
% holefilling.m
% Uwagi:
% do poprawnego dzia³ania funkcji konieczne jest zainstalowanie biblioteki
% WFDB Toolbox
% sygna³y i anotacje musz¹ byæ zgodne z formatem u¿ywanym w WFDB Toolbox
% odkomentowana zosta³a kolejnoœæ z najlepszymi wynikiami. Aby u¿yæ innych
% nale¿y je odkomentowaæ.
% Autor: Antonina Pater, 2018
% Ostatnia modyfikacja: 1.01.19

tolerancja = 0.2; % 200ms
jqrs_ann = run_qrsdet_by_seg(ECG_RawData,HRVparams);

% wczytanie sygna³u i anotacji
% [signal,Fs,tm] = rdsamp(name);
% ann1=rdann(name,'atr');
gauss = gauss_r(signal, Fs, 6, 20);
% load(strcat(name,'gauss.mat')); 
% gqrs(name,[],[],[],[],'gqrs')
% gann = rdann(name,'gqrs'); 
vann = waveletdet(signal, Fs);
% load(strcat(name,'wavelet.mat'));
% wqrs(name);
% wann=rdann(name,'wqrs'); 
% sqrs(name);
% sann=rdann(name,'qrs');

% uzupe³nienie anotacji w kolejnoœci: waveletedet, gauss_r, gqrs, wqrs,
% sqrs
newann = jqrs_ann';
[newann, ~] = holefilling(newann, gauss, Fs, tolerancja, signal);
[newann, ~] = holefilling(newann, vann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, gann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, wann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, sann, Fs, tolerancja, signal);


% % uzupe³nienie anotacji w kolejnoœci: waveletedet, gauss_r, wqrs,
% % sqrs
% newann = vann;
% [newann, ~] = holefilling(newann, gauss, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, wann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, sann, Fs, tolerancja, signal);
% 
% % uzupe³nienie anotacji w kolejnoœci: waveletedet, gauss_r
% newann = vann;
% [newann, ~] = holefilling(newann, gauss, Fs, tolerancja, signal);
% 
% % uzupe³nienie anotacji w kolejnoœci: wqrs, waveletedet, gauss_r, gqrs, 
% % sqrs
% newann = wann;
% [newann, ~] = holefilling(newann, vann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, gauss, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, gann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, sann, Fs, tolerancja, signal);
% 
% % uzupe³nienie anotacji w kolejnoœci: sqrs, waveletedet, gauss_r, gqrs, 
% % wqrs
% newann = sann;
% [newann, ~] = holefilling(newann, vann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, gann, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, gauss, Fs, tolerancja, signal);
% [newann, ~] = holefilling(newann, wann, Fs, tolerancja, signal);
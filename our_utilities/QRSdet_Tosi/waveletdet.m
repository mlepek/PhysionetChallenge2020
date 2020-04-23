% function [vann] = waveletdet(name, signal, tm)
function vann = waveletdet(signal, Fs)
% Przeznaczenie:
% funkcja s�u�y do detekcji zespo��w QRS za pomoc� wsp�czynnik�w
% zwracannych przez transformat� falkow�
% argumenty wej�ciowe:
% name - nazwa sygna�u
% Funkcja zwraca:
% vann - indeksy anotacji
% U�ywane funkcje:
% QRSdetect_new.m
% uwagi:
% do poprawnego dzia�ania funkcji konieczne jest zainstalowanie bibliotek
% WFDB Toolbox i Wavelet Toolbox
% sygna�y i anotacje musz� by� zgodne z formatem u�ywanym w WFDB Toolbox
% Autor: Antonina Pater, 2018
% Ostatnia modyfikacja: 1.01.19


nloops = 3; %LICZBA PETLI
tolerancja = 0.2;

% [signal,Fs,tm] = rdsamp(name);
%na razie tylko 1 kana��
signal = signal - movmean(signal,60);
signal = signal*0.7;
% ann1=rdann(name,'atr'); %do sprawdzania
k = 32 - mod(length(signal(:,1)),32);
signal = [signal(:,1); zeros(k,1)]; %swt chce �eby liczba pr�bek sygna�u by�a podzielna przez 2^5
% tm = [tm; zeros(k,1)];

[~,swd] = swt(signal,5,'haar'); %Haar wavelet transform
h = sqrt((abs((swd(4,:)).*(swd(5,:)))))';
L = 0.3 *max(h); %treshold
ann = find(h>= L); %indeksy nowowyktytych QRS�w

ann = [1 ; ann(diff(ann)>=0.2*Fs) ; length(h)];
RR = diff(ann); %interwa�y
newerann = ann;
for s =1:nloops
    [newerann, L]  = QRSdetect_new(RR, h, L, Fs, newerann, tolerancja, signal);
    RR = diff(newerann); %interwa�y
end

vann = newerann;
vann(vann==1) = [];
vann(vann >= length(signal)) = [];
% save(strcat(name,'wavelet.mat'),'vann');

end
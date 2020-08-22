% function hann = hierarchiczny_12(name, signal,Fs)

function hann = hierarchiczny_12(signal,Fs)
[no_samp, no_chann] =  size(signal);
window = round(0.1*Fs); % 100 ms
minnovotes = no_chann/2 + 1; % z ilu kanalow jest potrzebna minimalna liczba detekcji

% cell array na anotacje z wszystkich kanalow
ann_cell{1,no_chann} = [];

for channel = 1:no_chann
    % disp([num2str(channel), '/', num2str(no_chann)])
    ann_cell{1,channel} = hierarchiczny(signal(:,channel), Fs);
    ann_cell{1,channel}(ann_cell{1,channel}>length(signal)...
        | ann_cell{1,channel}<=window) = [];
end

% tablica do porównania anotacji
ann_comp = zeros(no_samp, no_chann);

for i=1:window*2
    for channel = 1:no_chann
        ann_comp((ann_cell{1,channel} - window + i), channel) = 1;
    end
end

%zliczanie "glosow"
ann_comp = ann_comp(1:no_samp, :);
sum_ann = sum(ann_comp,2);
new_ann = find(sum_ann >= minnovotes);
annDiff = diff(new_ann);
annDiffAbove0 = [1 ; find(annDiff>1)];

annNew = [];
for j=1:length(annDiffAbove0)-1
    annNew = [annNew ; round(mean(new_ann(annDiffAbove0(j):annDiffAbove0(j+1))))];
end
if any(ann_cell{1,1}) 
    hann = [annNew ; ann_cell{1,1}(end)];
else
   hann = annNew;
end

% save(strcat(name,'hann12.mat'),'hann');
end
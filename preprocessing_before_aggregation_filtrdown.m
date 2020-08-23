function data = preprocessing_before_aggregation_filtrdown(data, fs)

fs_fixed = 100; %docelowe sample frequency
%max_length = round(10*fs_fixed); %ustawienie maksymalnej dlugosci sygnalu: 10 s

%DOWNSAMPLING
    data= resample(data',fs_fixed,fs); %resamoling do 100Hz
    data = data';
    
    %DATA FILTERING
    ECG12filt = [];
    % median filter to remove bw
    for ii=1:12
        ECG12filt(ii,:) = movmean(data(ii,:)', 5);
    end
    
    ECG12filt2 = ECG12filt;
    for iii=1:size(ECG12filt,1)
        ECG12filt2(iii,:) = ECG12filt(iii,:) - movmedian(ECG12filt(iii,:), 50);
    end
    data = ECG12filt2;
%     
%     %PRZYCINANIE DLUZSZYCH SYGNALOW
%     matrixToExport_withzeros = zeros(12,max_length);
% %     if(size(data,2)>size(matrixToExport_withzeros,2)) %usinanie dluzszych sygnalow
% %         matrixToExport_withzeros(1:12,1:size(matrixToExport_withzeros,2)) = data(:,1:size(matrixToExport_withzeros,2));
% %     else
% %         matrixToExport_withzeros(1:12,1:size(data,2)) = data;
% %     end
% %     %             all_signals_matrix = [all_signals_matrix reshape(matrixToExport_withzeros',[],1)];
%     
%     data_new = [];
% 
%     if(size(data,2)>size(matrixToExport_withzeros,2)) %usinanie dluzszych sygnalow
%         matrixToExport_withzeros(1:12,1:size(matrixToExport_withzeros,2)) = data(:,1:size(matrixToExport_withzeros,2));
%     else
%         %kopiowanie za krotkich sygnalow:
%         if(size(data,2)<size(matrixToExport_withzeros,2))
%             howManyTimes = floor(size(matrixToExport_withzeros,2)/size(data,2)); %jaka jest wielokrotnosc dlugosci sygnalu (zaokraglona w dol)
%             howManyToCopy_residuum = size(matrixToExport_withzeros,2)-size(data,2)*howManyTimes; %reszta brakujacych probek
%             for j=1:howManyTimes
%                 data_new = [data_new data];
%             end
%            
%             %jesli zostaly jakies resztki (po dodaniu wielokrotnosci)
%             if(howManyToCopy_residuum>0)
%                 data_new = [data_new data(:,howManyToCopy_residuum)];
%             end
%             data = data_new;   
%         end
%        
%         matrixToExport_withzeros(1:12,1:size(data,2)) = data;
%     end
% 
% 
%     all_signals_matrix_output = reshape(matrixToExport_withzeros',[],1);
%     
%     
%     
%     
%     
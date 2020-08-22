function [RRMatrixToExport_withzeros,ann] = getRR(max_length,fs_fixed, data)

  addpath('./QRSdet_Tosi');
  
%   matrixToExport_withzeros = zeros(12,max_length);
%     coherentMatrixToExport = zeros(12, 1.2*fs_fixed); 
    RRMatrixToExport_withzeros = zeros(1, max_length+1.2*fs_fixed);
    %PRZYCINANIE DLUZSZYCH RRów
    RR_new = [];
    ann = hierarchiczny_12(data',fs_fixed);
    RR = diff(ann)/fs_fixed;
    if ~isempty(RR)
        if(length(RR)>length(RRMatrixToExport_withzeros)) %ucinanie dluzszych sygnalow
%             datatmp = data(:,1:size(matrixToExport_withzeros,2));
            RR = RR(1:size(RRMatrixToExport_withzeros,2));
            
        else
            %kopiowanie za krotkich sygnalow:
            if(length(RR)<length(RRMatrixToExport_withzeros))
                howManyTimes = floor(size(RRMatrixToExport_withzeros,2)/length(RR)); %jaka jest wielokrotnosc dlugosci sygnalu (zaokraglona w dol)
                howManyToCopy_residuum = size(RRMatrixToExport_withzeros,2)-length(RR)*howManyTimes; %reszta brakujacych probek
                
                for j=1:howManyTimes
                    RR_new = [RR_new; RR];
                end
                
                %jesli zostaly jakies resztki (po dodaniu wielokrotnosci)
                if(howManyToCopy_residuum>0)
                    RR_new = [RR_new; RR(1:howManyToCopy_residuum)];
                end
                RR = RR_new;
            end
%             ECG_header.nsamp = length(data);
            
        end
        
        RRMatrixToExport_withzeros(:, 1:size(RR,1)) = RR*1000; % zeby bylo w ms
    end
    
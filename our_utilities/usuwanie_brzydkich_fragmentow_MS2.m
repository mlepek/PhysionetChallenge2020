function [ECG12filt_cleaned, ECG12filt_binary] = usuwanie_brzydkich_fragmentow_MS2(ECG12filt, fs)
    
    len_sig = length(ECG12filt(1,:));
    
    ECG12filt_filtered = zeros(size(ECG12filt)); 
    ECG12filt_cleaned = ECG12filt;
    ECG12filt_binary = ones(size(ECG12filt));
    
    f1 = 5;
    f2 = 15;
    order = 2;
  
    [B,A] = butter(order,[f1/(fs/2) f2/(fs/2)]);
 
    for ch = 1:12
        signal = ECG12filt(ch,:);
%         signal=filter(B,A,signal);
%         ECG12filt_filtered(ch,:)=filter(B,A,signal1);  
%         signal = ECG12filt_filtered(ch,:);
        
        STD_overall = std(signal);
        kurtosis_overall = kurtosis(signal);
        
        STD_mov=movstd(signal,2*fs);

        ECG12filt_binary(ch, STD_mov>1.5*STD_overall) = 0;
        ECG12filt_cleaned(ch,STD_mov>1.5*STD_overall) = 0;
        
        
%         for i = fs:len_sig-fs % co 2 sekundy czyli 1000 próbek(?) 
%             
%             part_of_signal = signal(:,i-fs+1:i+fs);
%             kurt= kurtosis(part_of_signal);
%             dev = std(part_of_signal);
%             
%             if kurt < 3.6 || dev > 1.5*STD_overall
%                 
%                 ECG12filt_binary(ch, i-fs+1:i+fs) = 0;
%                 ECG12filt_cleaned(ch,i-fs+1:i+fs) = 0;
%                 
%             end
%         end
    end 
end

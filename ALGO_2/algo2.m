function features = algo2(data, header_data)
    % Input: sygna³ i nag³ówek dla tego sygna³u
    % Output: wektor (macierz) z uœrednionym beatem
    % 
    % Body:
    % 
    % a) Usuwanie brzydkich fragmentów.
    % b) Filtracja
    % c) Wyznaczanie QRS
    % d) Obliczanie uœrednionego beatu
    % e) Ewentualne dodatkowe operacje przed przekazaniem sygna³u na wyjœcie.
    
    %% ------------- BEGIN CODE --------------

    % read number of leads, sample frequency and gain from the header.	
	[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

	HRVparams.Fs=Fs;
    HRVparams.PeakDetect.windows = floor(Total_time-1);
    HRVparams.windowlength = floor(Total_time);

    for i =1:num_leads
            Lead12wGain(i,:) = data(i,:)* gain(i);
    end

    % median filter to remove bw
    for i=1:num_leads
        ECG12filt(i,:) = movmedian(Lead12wGain(i,:)', 5);
    end

    [ECG12filt, ~] = usuwanie_brzydkich_fragmentow_MS2(ECG12filt, HRVparams.Fs);

    ECG12filt2 = ECG12filt;
    for i=1:size(ECG12filt,1)
        ECG12filt2(i,:) = ECG12filt(i,:) - movmedian(ECG12filt(i,:), HRVparams.Fs/8);
    end

    XYZLeads = Kors_git(ECG12filt2);
    VecMag_toQRS = vecnorm(XYZLeads');

    % convert 12Leads to XYZ leads using Kors transformation
    XYZLeads = Kors_git(ECG12filt);

    VecMag = vecnorm(XYZLeads');

    % Convert ECG waveform in rr intervals
    [t, rr, jqrs_ann, ~ , ~] = ConvertRawDataToRRIntervals(VecMag_toQRS, HRVparams, recording, ECG12filt');

    % Find fiducial points using ECGKit
    ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
    wavedet_config.setup.wavedet.QRS_detection_only = 0;
    [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
     
%% 
    [features,Fid_pts_Median] = Time_coherent_code_github_algo2(ECG12filt2', Fid_pts, Fs);
        
end

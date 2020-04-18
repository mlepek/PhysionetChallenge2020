function features = get_12ECG_features(data, header_data)

       % addfunction path needed
        addpath(genpath('Tools/'))
        load('HRVparams_12ECG','HRVparams')

        

	% read number of leads, sample frequency and gain from the header.	

	[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

	HRVparams.Fs=Fs;
        HRVparams.PeakDetect.windows = floor(Total_time-1);
        HRVparams.windowlength = floor(Total_time);
% tic
	try

                for i =1:num_leads
                        Lead12wGain(i,:) = data(i,:)* gain(i);
                end


                % median filter to remove bw
                for i=1:num_leads
                        ECG12filt(i,:) = medianfilter(Lead12wGain(i,:)', Fs);
                end

                % convert 12Leads to XYZ leads using Kors transformation
                XYZLeads = Kors_git(ECG12filt);

                VecMag = vecnorm(XYZLeads');


                % Convert ECG waveform in rr intervals
                [t, rr, jqrs_ann, ~ , ~] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
%                 sqi = [tSQI', SQIvalue'];

                % Find fiducial points using ECGKit
                ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
                wavedet_config.setup.wavedet.QRS_detection_only = 0;
                [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);

                [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);
   
                
                GEH_features = GEH_analysis_git(XYZ_Median,Fid_pts_Median,Fs);

                features(1)=age;
                features(2)=sex;
                features(3:24)=GEH_features;
                

	catch
		features = NaN(1,24);
    end
% ttt_default = toc   
    %% Our own excellent features here
    
     addpath('getting_features/');
% tic    
    try
        % features_normal_muter
        
        % comment: prosba o ubranie w funkcje
        features_muter_all = get_function_Muter(Fid_pts, ECG12filt, Fs,header_data);

        % taking all Kasia Muter's features together        
        features(25:33) = features_muter_all;
    catch
        features(25:33) = NaN(1,9);
    end
% ttt_muter = toc    
    %
% tic
    try
        [TFRPrimaV1V2V3, RSmax, QRSwidth, SyRyMedianV1, QRpercmed, QRStimemin] = ...
                                    get_function_Kokosinska(ECG12filt,XYZLeads,header_data);                         
        features(33:38) = [TFRPrimaV1V2V3, RSmax, QRSwidth, SyRyMedianV1, QRpercmed, QRStimemin];
    catch
        features(33:38) = NaN(1,6);
    end
    %}
% ttt_koko = toc    
    %
    
% tic
    try
        
        [PACfeatures, PVCfeatures] = get_function_Pater(XYZLeads, ECG12filt(2,:), Fs,header_data);  
        
        PACcell  = struct2cell(PACfeatures);
        empties = cellfun('isempty',PACcell);
        PACcell(empties) = {NaN};
        PACvec = cat(2,PACcell{:});
        PVCcell = struct2cell(PVCfeatures);
        empties = cellfun('isempty',PVCcell);
        PVCcell(empties) = {NaN};
        PVCvec = cat(2,PVCcell{:});
        features(39:51) = [PACvec, PVCvec];
    catch
        features(39:51) = NaN(1,13);
    end
    %}
% ttt_pater = toc


% tic
    try
        features_wiszniewski = get_function_wiszniewski(XYZ_Median, Fs, Fid_pts, Fid_pts_Median, rr);                      
        features(52:56) = features_wiszniewski;
    catch
        features(52:56) = NaN(1,5);
    end
% ttt_przemek = toc


end

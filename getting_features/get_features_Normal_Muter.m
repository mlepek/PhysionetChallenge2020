function features = get_features_Normal_Muter(data, header_data)

% Input:
% - data
% - header_data
% 
% Output features:
% - œrednia wartoœæ heart rate (*Regular rhythm at a rate of 60-100 bpm)
% - procent za³amków P wzglêdem za³amków R
% - znak P w odpr. I (1 - jeœli >= 0; 0 - jeœli <0)
% - znak P w odpr. II (1 - jeœli >= 0; 0 - jeœli <0)
% - znak P w odpr. aVR (1 - jeœli >= 0; 0 - jeœli <0)
% - odch. std. odcinków PR

       % addfunction path needed
        addpath(genpath('Tools/'))
        load('HRVparams_12ECG','HRVparams')

	% read number of leads, sample frequency and gain from the header.	

	[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

	HRVparams.Fs=Fs;
        HRVparams.PeakDetect.windows = floor(Total_time-1);
        HRVparams.windowlength = floor(Total_time);

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
                [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
%                 sqi = [tSQI', SQIvalue'];

                % Find fiducial points using ECGKit
                ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
                wavedet_config.setup.wavedet.QRS_detection_only = 0;
                [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);

%                 [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);

                % Each QRS complex is preceded by a normal P wave
                
                detected_P_waves = length(Fid_pts.P(~isnan(Fid_pts.P)));
                detected_R_waves = length(Fid_pts.R(~isnan(Fid_pts.R)));
                percent_of_P = (detected_P_waves/detected_R_waves)*100;
                
                % Heart Rate, fs = 500 Hz
                
                RR_interval = [];
                
                for i = 1:detected_R_waves-1
                    if (isnan(Fid_pts.R(i))==true || isnan(Fid_pts.R(i+1))==true)
                        continue;
                    else
                        RR_temp = Fid_pts.R(i+1) - Fid_pts.R(i);
                    end
                    RR_interval = [RR_interval, RR_temp];
                end
                
                RR_mean = mean(RR_interval);                
                HR = 60/(RR_mean/Fs);
                
                % Normal P wave axis: 
                % - P waves should be upright in leads I and II, 
                % - inverted in aVR.
                
                P_waves = Fid_pts.P(:,~isnan(Fid_pts.P));
                
                mean_P_L1 = 0;
                mean_P_L2 = 0;
                mean_P_aVR = 0;
                
                for i=1:detected_P_waves
                    
                    xL1 = ECG12filt(1,P_waves(i));
                    xL2 = ECG12filt(2,P_waves(i));
                    xaVR = ECG12filt(4,P_waves(i));
                    
                    mean_P_L1 = mean_P_L1 + xL1;
                    mean_P_L2 = mean_P_L2 + xL2;
                    mean_P_aVR = mean_P_aVR + xaVR;
                    
                end
                mean_P_L1 = mean_P_L1/detected_P_waves;
                mean_P_L2 = mean_P_L2/detected_P_waves;
                mean_P_aVR = mean_P_aVR/detected_P_waves;
                                         
                % The PR interval and P wave morphology in the given lead remain constant.
                
                PR_interval = [];
                
                for i = 1:detected_R_waves
                    if (isnan(Fid_pts.P(i))==true || isnan(Fid_pts.R(i))==true)
                        continue;
                    else
                        PR_temp = Fid_pts.R(i) - Fid_pts.P(i);
                    end
                    
                    PR_interval = [PR_interval, PR_temp];
                end
                
                std_PR = std(PR_interval);
                
                % Features vector
                
                features(1) = HR;       % sredni rytm serca
                features(2) = percent_of_P;     % procent za³. P wzglêdem R
                if mean_P_L1>=0         % znak P w odpr. I
                    features(3) = 1;
                else
                    features(3) = 0;
                end
                if mean_P_L2>=0         % znak P w odpr. II
                    features(4) = 1;
                else
                    features(4) = 0;
                end
                if mean_P_aVR>=0        % znak P w odpr. aVR
                    features(5) = 1;
                else
                    features(5) = 0;
                end
                features(6) = std_PR;   % odch. std. odcinka PR


	catch
		features = NaN(1,6);
	end

end


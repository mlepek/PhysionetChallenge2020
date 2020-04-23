function features = get_features_STD_Muter(data, header_data)

% STD - Obni¿enie odcinka ST,
%       Horizontal or downsloping ST depression  0.5mm
% 
% Input:
% - data
% - header_data
% 
% Output features:
% Max. amplituda miedzy J (QRSoffset) a J+80 ms w odprowadzeniach:
% - II, III, aVF
% - V1-V4
% - I, aVL, V5 i V6
% (Obliczamy cechê na uœrednionych beatach w trzech poszczególnych grupach
% odprowadzeñ)

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
                sqi = [tSQI', SQIvalue'];

                % Find fiducial points using ECGKit
                ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
                wavedet_config.setup.wavedet.QRS_detection_only = 0;
                [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);

%                 [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);
                
                
                % Amplituda miêdzy J (QRSoffset) a J+80 ms
                
                % II, III, aVF
                
                K2 = [ECG12filt(2,:);ECG12filt(2,:);ECG12filt(2,:)]';
                [K2_Median,Fid_pts_K2_Median] = Time_coherent_code_github(K2,Fid_pts,Fs);
                amplK2 = K2_Median(Fid_pts_K2_Median.QRSoff,1) - K2_Median(Fid_pts_K2_Median.QRSoff+40,1);
                
                K3 = [ECG12filt(3,:);ECG12filt(3,:);ECG12filt(3,:)]';
                [K3_Median,Fid_pts_K3_Median] = Time_coherent_code_github(K3,Fid_pts,Fs);
                amplK3 = K3_Median(Fid_pts_K3_Median.QRSoff,2) - K3_Median(Fid_pts_K3_Median.QRSoff+40,1);
                                
                aVF = [ECG12filt(6,:);ECG12filt(6,:);ECG12filt(6,:)]';
                [aVF_Median,Fid_pts_aVF_Median] = Time_coherent_code_github(aVF,Fid_pts,Fs);
                amplaVF = aVF_Median(Fid_pts_aVF_Median.QRSoff,3) - aVF_Median(Fid_pts_aVF_Median.QRSoff+40,1);
                
                ampl_K2K3aVF_max = max([amplK2, amplK3, amplaVF]);
                
                % V1-V4
                
                V1 = [ECG12filt(7,:);ECG12filt(7,:);ECG12filt(7,:)]';
                [V1_Median,Fid_pts_V1_Median] = Time_coherent_code_github(V1,Fid_pts,Fs);
                amplV1 = V1_Median(Fid_pts_V1_Median.QRSoff,1) - V1_Median(Fid_pts_V1_Median.QRSoff+40,1);
                
                V2 = [ECG12filt(8,:);ECG12filt(8,:);ECG12filt(8,:)]';
                [V2_Median,Fid_pts_V2_Median] = Time_coherent_code_github(V2,Fid_pts,Fs);
                amplV2 = V2_Median(Fid_pts_V2_Median.QRSoff,1) - V2_Median(Fid_pts_V2_Median.QRSoff+40,1);
                
                V3 = [ECG12filt(9,:);ECG12filt(9,:);ECG12filt(9,:)]';
                [V3_Median,Fid_pts_V3_Median] = Time_coherent_code_github(V3,Fid_pts,Fs);
                amplV3 = V3_Median(Fid_pts_V3_Median.QRSoff,1) - V3_Median(Fid_pts_V3_Median.QRSoff+40,1);
                
                V4 = [ECG12filt(10,:);ECG12filt(10,:);ECG12filt(10,:)]';
                [V4_Median,Fid_pts_V4_Median] = Time_coherent_code_github(V4,Fid_pts,Fs);
                amplV4 = V4_Median(Fid_pts_V4_Median.QRSoff,1) - V4_Median(Fid_pts_V4_Median.QRSoff+40,1);
                            
                ampl_V1V2V3V4_max = max([amplV1, amplV2, amplV3, amplV4]);
                
                % I, aVL,  V5 i V6
                
                K1 = [ECG12filt(1,:);ECG12filt(1,:);ECG12filt(1,:)]';
                [K1_Median,Fid_pts_K1_Median] = Time_coherent_code_github(K1,Fid_pts,Fs);
                amplK1 = K1_Median(Fid_pts_K1_Median.QRSoff,1) - K1_Median(Fid_pts_K1_Median.QRSoff+40,1);
                
                aVL = [ECG12filt(5,:);ECG12filt(5,:);ECG12filt(5,:)]';
                [aVL_Median,Fid_pts_aVL_Median] = Time_coherent_code_github(aVL,Fid_pts,Fs);
                amplaVL = aVL_Median(Fid_pts_aVL_Median.QRSoff,3) - aVL_Median(Fid_pts_aVL_Median.QRSoff+40,1);
                
                V5 = [ECG12filt(11,:);ECG12filt(11,:);ECG12filt(11,:)]';
                [V5_Median,Fid_pts_V5_Median] = Time_coherent_code_github(V5,Fid_pts,Fs);
                amplV5 = V5_Median(Fid_pts_V5_Median.QRSoff,1) - V5_Median(Fid_pts_V5_Median.QRSoff+40,1);
                
                V6 = [ECG12filt(12,:);ECG12filt(12,:);ECG12filt(12,:)]';
                [V6_Median,Fid_pts_V6_Median] = Time_coherent_code_github(V6,Fid_pts,Fs);
                amplV6 = V6_Median(Fid_pts_V6_Median.QRSoff,1) - V6_Median(Fid_pts_V6_Median.QRSoff+40,1);
                        
                ampl_K1aVLV5V6_max = max([amplK1, amplaVL, amplV5, amplV6]);
                
                % Features vector
                features(1) = ampl_K2K3aVF_max;
                features(2) = ampl_V1V2V3V4_max;
                features(3) = ampl_K1aVLV5V6_max;

	catch
		features = NaN(1,3);
	end

end


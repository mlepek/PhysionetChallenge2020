function [baselinePR_mean, baselineT, min_baseline] = get_function_baseline(ECG12filt,header_data, lead_num)
addpath(genpath('Tools/'))
load('HRVparams_12ECG','HRVparams')
baselineT = 0;
[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

HRVparams.Fs=Fs;
HRVparams.PeakDetect.windows = floor(Total_time-1);
HRVparams.windowlength = floor(Total_time);

% Toff = [Fid_pts_Median(:).Toff];
% TuV=XYZ_Median(Toff);
% h=figure;
% plot(XYZ_Median)
% hold on;
% yline(TuV,'-.b')
% PR = [Fid_pts_Median(:).QRSoff];
% PRuV=XYZ_Median(PR);
% h=figure;
% plot(XYZ_Median)
% hold on;
% yline(PRuV,'-.b')

%Œredni beat w we wszystkich odprowadzeniach
% leadsall=[1 2 3 4 6 7 8 9 10 11 12]; %bez aVL,bo coœ jest nie tak :/
% for k=1:length(leadsall)
    i = lead_num;      
    ALL = [ECG12filt(i,:);ECG12filt(i,:);ECG12filt(i,:)]';
    VecMag = vecnorm(ALL');
    [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
%     sqi = [tSQI', SQIvalue'];
    ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
    wavedet_config.setup.wavedet.QRS_detection_only = 0;
    [Fid_pts_ALL,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
%     [XYZ_Median_ALL,Fid_pts_Median_ALL] = Time_coherent_code_github(ALL,Fid_pts_ALL,Fs);
   
    baselinePR_mean = mean(ECG12filt(i,Fid_pts_ALL.QRSon(~isnan(Fid_pts_ALL.QRSon)))); %srednia ze wszystkich QRSon w odprowadzeniu
    %Tend------------------------------------------------------------------
    % SAI QRST was calculated as the arithmetic sum of areas under the QRST curve on XYZ leads, with baseline defined as the voltage at the end of the T-wave. 
%     Toff = [Fid_pts_Median_ALL(:).Toff];
%     baselineT=XYZ_Median_ALL(Toff);
    
    baselineT_mean = mean(ECG12filt(i,Fid_pts_ALL.Toff(~isnan(Fid_pts_ALL.Toff)))); %srednia ze wszystkich QRSon w odprowadzeniu

%     h(k)=figure;
%     plot(XYZ_Median_ALL);
%     hold on;
%     plotbaselineT=yline(baselineT,'-.b');
    % Create a PNG filename.
%     pngFileName =[recording,'T_',num2str(i),'.png']
%     fullFileName = fullfile('C:\Users\dkokosinska\Desktop\Physionet Challenge\Training_WFDB\Spr\Output\Baseline_figures', pngFileName);
    
    % Then save it
%     export_fig(fullFileName);
    
    
    %PRsegment----------------------------------------------------------
%     PR = [Fid_pts_Median_ALL(:).QRSon];
%     baselinePR=XYZ_Median_ALL(PR);
%     figure(1);
% %     plot(XYZ_Median_ALL);
%     plot(ECG12filt(i,:));
%     hold on;
% %     plotbaselinePR=yline(baselinePR,'-.b');
%     plotbaselineT=yline(baselineT_mean,'-.m');
%     plotbaselinePR_mean = yline(baselinePR_mean,'-.g');
%     plot(Fid_pts_Median_ALL.QRSon,XYZ_Median_ALL(Fid_pts_Median_ALL.QRSon),'r.');
%     hold off
    % Create a PNG filename.
%     pngFileName =[recording,'PR_',num2str(i),'.png']
%   	fullFileName = fullfile('C:\Users\dkokosinska\Desktop\Physionet Challenge\Training_WFDB\Spr\Output\Baseline_figures', pngFileName);
    % Then save it
%     export_fig(fullFileName);
    %Jpoint
%     Jpoint = [Fid_pts_Median_ALL(:).QRSoff];
%     JpointuV=XYZ_Median_ALL(Jpoint);
    min_baseline = min([baselinePR_mean baselineT_mean]);
% end
end






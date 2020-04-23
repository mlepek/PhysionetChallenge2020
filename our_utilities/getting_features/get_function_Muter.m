function features = get_function_Muter(Fid_pts, ECG12filt, Fs, header_data)
%  Input:,
%
% - Fid_pts - struct vector with the detected points locations in samples
%              from wavedet_3D_ECGKit function
% - ECG12filt - filtered ECG signals
% - Fs - sampling frequency
%
% Output features - Normal:
% - œrednia wartoœæ heart rate (*Regular rhythm at a rate of 60-100 bpm)
% - procent za³amków P wzglêdem za³amków R
% - znak P w odpr. I (1 - jeœli >= 0; 0 - jeœli <0)
% - znak P w odpr. II (1 - jeœli >= 0; 0 - jeœli <0)
% - znak P w odpr. aVR (1 - jeœli >= 0; 0 - jeœli <0)
% - odch. std. odcinków PR

% STD - Obni¿enie odcinka ST,
%       Horizontal or downsloping ST depression  0.5mm

% Output features - STD:
% Max. amplituda miedzy J (QRSoffset) a J+80 ms w odprowadzeniach:
% - II, III, aVF
% - V1-V4
% - I, aVL, V5 i V6
% (Obliczamy cechê na uœrednionych beatach w trzech poszczególnych grupach
% odprowadzeñ)


% -------------------------NORMAL------------------------

% Each QRS complex is preceded by a normal P wave
if ~isempty(Fid_pts)
    detected_P_waves = length(Fid_pts.P(~isnan(Fid_pts.P)));
    detected_R_waves = length(Fid_pts.R(~isnan(Fid_pts.R)));
    if(detected_R_waves>0)
        percent_of_P = (detected_P_waves/detected_R_waves)*100;
    elseif (detected_R_waves==0)
        percent_of_P = NaN;
        %         disp("percent_of_P = NaN: detected_R_waves = " + detected_R_waves + " detected_P_waves = " +detected_P_waves)
    else
        percent_of_P = 0;
    end
else
    percent_of_P = NaN;
    %     disp("percent_of_P = NaN: Fid_pts.R = " + Fid_pts.R + " Fid_pts.P = " +Fid_pts.P)
end
% Heart Rate, fs = 500 Hz

RR_interval = [];

RR_interval = diff(Fid_pts.R(~isnan(Fid_pts.R)));
if ~isempty(RR_interval)
    RR_mean = mean(RR_interval);
    HR = 60/(RR_mean/Fs);
else
    HR = NaN;
    %     disp("HR = NaN: RR_interval is empty.")
end

% Normal P wave axis:
% - P waves should be upright in leads I and II,
% - inverted in aVR.

P_waves = Fid_pts.P(:,~isnan(Fid_pts.P));

xL1 = ECG12filt(1,P_waves);
xL2 = ECG12filt(2,P_waves);
xaVR = ECG12filt(4,P_waves);

if ~isempty(P_waves)
    mean_P_L1 = mean(xL1);
    mean_P_L2 = mean(xL2);
    mean_P_aVR = mean(xaVR);
else
    %     disp("P_waves is empty.")
    mean_P_L1 = NaN;
    mean_P_L2 = NaN;
    mean_P_aVR = NaN;
end

% The PR interval and P wave morphology in the given lead remain constant.

PR_interval = Fid_pts.R(~isnan(Fid_pts.P) & ~isnan(Fid_pts.R)) - Fid_pts.P(~isnan(Fid_pts.P) & ~isnan(Fid_pts.R));
if (~isempty(PR_interval) && length(PR_interval) > 1)
    std_PR = std(PR_interval);
elseif length(PR_interval)==1
    std_PR = NaN;
    %     disp("Length of PR_interval is 1, std is NaN.")
else
    std_PR = NaN;
    %     disp("std_PR is NAN: Fid_pts.R = " + num2str(Fid_pts.R) + " Fid_pts.P = " +num2str(Fid_pts.P)+" PR_interval is empty.")
end

% -------------------------- STD ----------------------------

% Amplituda  J+80 ms - a baseline (doda³em oblicznie baseline - matsol
% 10042020)
% tic
% II, III, aVF

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,2);
    amplK2 = mean(ECG12filt(2,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplK2 = NaN;
end

% K2 = [ECG12filt(2,:);ECG12filt(2,:);ECG12filt(2,:)]';
% [K2_Median,Fid_pts_K2_Median] = Time_coherent_code_github(K2,Fid_pts,Fs);
% % amplK2 = K2_Median(Fid_pts_K2_Median.QRSoff,1) - K2_Median(Fid_pts_K2_Median.QRSoff+40,1);
% amplK2 = K2_Median(Fid_pts_K2_Median.QRSoff+40,1)-min_baseline;

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,3);
    amplK3 = mean(ECG12filt(3,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplK3 = NaN;
end

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,6);
    amplaVF = mean(ECG12filt(6,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplaVF = NaN;
end

ampl_K2K3aVF_max = max([amplK2, amplK3, amplaVF]);

% V1-V4

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,7);
    amplV1 = mean(ECG12filt(7,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplV1 = NaN;
end

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,8);
    amplV2 = mean(ECG12filt(8,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplV2 = NaN;
end

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,9);
    amplV3 = mean(ECG12filt(9,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplV3 = NaN;
end

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,10);
    amplV4 = mean(ECG12filt(10,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplV4 = NaN;
end

ampl_V1V2V3V4_max = max([amplV1, amplV2, amplV3, amplV4]);

% I, aVL,  V5 i V6

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,1);
    amplK1 = mean(ECG12filt(1,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplK1 = NaN;
end

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,5);
    amplaVL = mean(ECG12filt(5,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplaVL = NaN;
end

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,11);
    amplV5 = mean(ECG12filt(11,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplV5 = NaN;
end

try
    [~, ~,min_baseline] = get_function_baseline(ECG12filt,header_data,12);
    amplV6 = mean(ECG12filt(12,Fid_pts.QRSoff(~isnan(Fid_pts.QRSoff))+40,1)-min_baseline);
catch
    amplV6 = NaN;
end

ampl_K1aVLV5V6_max = max([amplK1, amplaVL, amplV5, amplV6]);
% toc
% ----------------------- FEATURES VECTOR ----------------------------

% Features vector - Normal

features(1) = HR;       % sredni rytm serca
features(2) = percent_of_P;     % procent za³. P wzglêdem R
if mean_P_L1>=0         % znak P w odpr. I
    features(3) = 1;
elseif mean_P_L1<0
    features(3) = 0;
else
    features(3) = NaN;
end

if (mean_P_L2 >= 0)         % znak P w odpr. II
    features(4) = 1;
elseif (mean_P_L2 < 0)
    features(4) = 0;
else
    features(4) = NaN;
end

if (mean_P_aVR >= 0)        % znak P w odpr. aVR
    features(5) = 1;
elseif (mean_P_aVR < 0)
    features(5) = 0;
else
    features(5) = NaN;
end

features(6) = std_PR;   % odch. std. odcinka PR

% Features vector - STD

if ~isnan(ampl_K2K3aVF_max)
    features(7) = ampl_K2K3aVF_max;
else
    features(7) = NaN;
%     disp("ampl_K2K3aVF_max is NAN")
end

if ~isnan(ampl_V1V2V3V4_max)
    features(8) = ampl_V1V2V3V4_max;
else
    features(8) = NaN;
%     disp("ampl_V1V2V3V4_max is NAN")
end

if ~isnan(ampl_K1aVLV5V6_max)
    features(9) = ampl_K1aVLV5V6_max;
else
    features(9) = NaN;
%     disp("ampl_K1aVLV5V6_max is NAN")
end

end
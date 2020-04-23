function features = get_function_wiszniewski(XYZ_Median, Fs, Fid_pts, Fid_pts_Median, RR)
% %Input:
%     XYZ_Median - averaged ECG signal from X, Y, Z leads
%     Fs - sampling frequency
%     Fid_pts - struct vector with the detected points locations in samples
%               from wavedet_3D_ECGKit function
%     Fid_pts_Median - median Fid points
%     RR - rr intervals vector from ConvertRawDataToRRIntervals function
% 
% %Output:
%     Vector containing features:
%       - P to R ratio
%       - mean area under fragment XYZ_Median(QRSonset-40 samples, QRSonset)
%         and mean of the absolute slopes of the fitted lines to this fragment 
%         for each X,Y,Z leads
%       - STD(RR)
%       - RMSSD(RR)
%       - mean PR interval

%% AF
%ta cecha pojawia siê ju¿ u Kasi
% P_R_ratio = nnz(~isnan(Fid_pts.P)) / nnz(~isnan(Fid_pts.R));

num_prev_samp = 40;
samp_range = Fid_pts_Median.QRSon-num_prev_samp:Fid_pts_Median.QRSon-1;
frag_bef_QRS = XYZ_Median(samp_range, :);
area_bef_QRS = sum(abs(frag_bef_QRS - mean(frag_bef_QRS))) / Fs;
c_x = polyfit(samp_range', frag_bef_QRS(:, 1), 1);
c_y = polyfit(samp_range', frag_bef_QRS(:, 2), 1);
c_z = polyfit(samp_range', frag_bef_QRS(:, 3), 1);

%zmieni³em mean na max
mean_slope = max(abs([c_x(1) c_y(1) c_z(1)])) * Fs;
% RR = diff(Fid_pts.R) / Fs;    % unnecessary
std_RR = std(RR);
rmssd_RR = rms(diff(RR));

%% I-AVB
mean_PR = nanmean(diff([Fid_pts.P; Fid_pts.QRSon]) / Fs);
features = [mean(area_bef_QRS) mean_slope std_RR rmssd_RR mean_PR];
end


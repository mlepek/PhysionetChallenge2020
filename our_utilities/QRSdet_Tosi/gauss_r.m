function [peak_candidates] = gauss_r(ecg, fs, lowfreq, highfreq)
%RPEAKDETECT
% 
%   Input:
%   ecg - raw ECG signal
%   fs - sampling frequency
%   lowfreq -  low frequency of the band pass filter, default: 6 Hz
%   highfreq - high frequency of the band pass filter, default: 20 Hz
% 
%   Output:
%   peakI - indices of R peaks (shifted)
%   
%   Example of usage:
%   [peakI] = rpeakdetect_pub(ecg, 300, 6, 20);
% 
%   Source:
%   ECG heart beat detection based on:
%   An Efficient R-peak Detection Based on New Nonlinear 
%   Transformation and First-Order Gaussian Differentiator
%   P. Kathirvel, M. Sabarimalai Manikandan
%   Cardiovascular Engineering and Technology, 2011, p. 408�425
% 	link: http://link.springer.com/article/10.1007/s13239-011-0065-3/fulltext.html
%   Score (in paper): sensitivity of 99.94%, and a positive predictivity 
%   of 99.96% - on MIT_BIH Arhhythmia Database.
%   
%   by Ania & Jacek
%   Last modification: 06.04.17

%% 1. Bandpass Filtering and Differentiation

    % 15-th FIR filter using the least squares approach
    n1 = 15;
    f = 0:fs/2-1;
    f = f/(fs/2);
    a = [zeros(1,lowfreq) ones(1,highfreq-lowfreq) zeros(1,(fs/2)-highfreq)];
    b1 = firls(n1,f,a);
    ecg_band = filtfilt(b1,1,ecg);
    
    % differentiation
    ecg_diff = diff(ecg_band);
    
%% 2.New Nonlinear Transformation

    % squaring
    ecg_s = ecg_diff.^2;
    
    % thresholding and normalization
    thresholds = [];

    window = fs*5;
    for i = 1:length(ecg_s)/window
		d = ecg_s(1 +(i-1)*window:(i)*window);
		thresholds = [thresholds 0.5*std(d)];
    end
	threshold = median(thresholds);

	ecg_s(ecg_s < threshold) = 0;
    ecg_s = ecg_s/max(abs(ecg_s));

    ecg_ss = ecg_s.^2;
    
    % Shannon energy
    SE = - ecg_ss.*(log(ecg_ss));
    SE(~isfinite(SE)) = 0;
    
    % filtering Shannon energy
    n2 = 15;
    b2 = firls(n2,f,a);
    SE_filt = filtfilt(b2,1,SE);
    
%% 3. New Peak-Finding Technique

    % gaussian filter design 
    sigma = floor(fs/8)-1;
    M = ceil(fs*2.5);
    x = linspace(1, M, M);
    gaussFilter = exp(-(x-M/2).^2 / (2 * sigma ^ 2));
    
    FOGD = diff(gaussFilter);
    
    ZC = conv(SE_filt, FOGD, 'same');
    
	ZCs = find((ZC(1:end-1) > 0) & (ZC(2:end) < 0));
	peak_candidates = ZCs;
    
    % sorry, poor removal of bad peak candidates 
    slope = ZC(ZCs+1) - ZC(ZCs);
    %peak_candidates(slope>-10e-15) = []; % 10e-15 - empirical parameter 
    peak_candidates(slope>mean(slope)/10e10) = []; % mean(slope)/10e10 - also empirical parameter, but less vulnerable to diffrent signals, as we think so 
    peak_candidates(diff(peak_candidates)<0.25*fs) = [];
    %peak_candidates(diff(peak_candidates)>1.5*fs) = [];
 
%% 4. Finding Location of Real R-Peaks
% This part wasn't working well for our signals and we decided it's unnecessary.

%     K = 30;
%     zpd = zeros(1, K/2);
%     
%     if peak_candidates(1) < K/2
%         peak_candidates = peak_candidates + K/2;
%         ecg_filt = [zpd ecg];
%     else
%         ecg_filt = [ecg zpd];
%     end
%     
%     for m = 1:length(peak_candidates)
%         SegStartPosition = peak_candidates -  K/2;
%         SegEndPosition = peak_candidates +  K/2;
%         ECGSeg = ecg_filt(SegStartPosition:SegEndPosition);
%         [~, max_idx] = max(abs(ECGSeg));
%         peak_candidates = max_idx + SegStartPosition;
%     end
    
end


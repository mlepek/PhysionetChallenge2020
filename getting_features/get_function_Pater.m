% INPUT
% XYZLeads
% sig2 - ecg from lead II
% Fs - sampling frequency

% OUTPUT
% PACfeatures:
% - p_len_std - std(d³ugoœci za³amka P)
% - p_len_max_min - Max-min(d³ugoœci zalamka P)
% - p_ampl_std - std(amplitudy za³amka P)
% - p_ampl_max_min - Max-min(amplitudy zalamka P)
% - p_area_std - std(pole powierzchni pod zalamkiem P)
% - maxRRi - maksymalna ró¿nica miêdzy dwoma s¹siednimi interwa³ami
% RR max(RR_i+1 - RR_i)
% PVCfeatures:
% - qrs_en_std - std(z energii zespo³ów QRS)
% - qrs_en_max_min - Max-min(z energii zespo³ów QRS)
% - qrs_len_std - std (d³ugoœæ zespo³u QRS)
% - qrs_len_max_min - Max-min(d³ugoœæ zespo³u QRS)
% - t_ampl_anyneg - Wartoœæ amplitudy T przynajmniej raz by³a ujemna
% w odprowadzeniu 1. (1-tak, 0-nie)
% - qrsoff_anyneg - Ró¿nica wartoœci pomiêdzy punktami QRSoffset
% (punkt J) - baseline przynajmniej raz by³a ujemna (lub mniejsza ni¿
% jak¹œ ma³¹ ujemn¹ wartoœæ)   (1-tak, 0-nie)
% - anypause - Warunek, ¿e
% RR_i+1>RR_i & (RR_i-1 >= 0.90*(RR_i+1+RR_i)/2 ||  RR_i-1<= 1.10*(RR_i+1+RR_i)/2)
% spe³niony przynajmniej 1 raz.


function [PACfeatures, PVCfeatures] = get_function_Pater( XYZLeads, sig2, Fs,header_data)

% fiducial points
sigX = XYZLeads(:,1)'; sigY = XYZLeads(:,2)'; sigZ = XYZLeads(:,3)';
try
    Fid_pts_2chann = fiducial_points(sig2, header_data); % for II lead
catch
    Fid_pts_2chann = [];
end
try
    Fid_pts_X = fiducial_points(sigX, header_data); % for X lead
catch
    Fid_pts_X = [];
end
try
    Fid_pts_Y = fiducial_points(sigY, header_data); % for Y lead
catch
    Fid_pts_Y = [];
end
try
    Fid_pts_Z = fiducial_points(sigZ, header_data); % for Z lead
catch
    Fid_pts_Z = [];
    
end

%PAC
if isempty(Fid_pts_2chann)
    
    PACfeatures.p_len_std = NaN;
    PACfeatures.p_len_max_min = NaN;
    PACfeatures.p_ampl_std = NaN;
    PACfeatures.p_ampl_max_min = NaN;
    PACfeatures.p_area_std = NaN;
    
else
    p_len = Fid_pts_2chann.Poff - Fid_pts_2chann.Pon;
    PACfeatures.p_len_std = std(p_len, 'omitnan')/Fs;
    %brakowalo nawiasow w (max-min)/fs - matsol 09042020
    PACfeatures.p_len_max_min = (max(p_len)-min(p_len))/Fs;
    
    p_ampl = sig2(Fid_pts_2chann.P(~isnan(Fid_pts_2chann.P)));
    PACfeatures.p_ampl_std = std(p_ampl, 'omitnan');
    PACfeatures.p_ampl_max_min = (max(p_ampl)-min(p_ampl))/Fs;
    
    Fid_pon_omitnan = Fid_pts_2chann.Pon(~(isnan(Fid_pts_2chann.Pon)|isnan(Fid_pts_2chann.Poff)));
    Fid_poff_omitnan = Fid_pts_2chann.Poff(~(isnan(Fid_pts_2chann.Poff)|isnan(Fid_pts_2chann.Pon)));
    p_area = [];
    for i = 1:length(Fid_pon_omitnan)
        p_frag = sig2(Fid_pon_omitnan(i):Fid_poff_omitnan(i));
        p_area = [p_area; trapz(p_frag)];
    end
    PACfeatures.p_area_std = std(p_area)/Fs;
end

%PVC
RR_interval = [];
qrs_enX = []; qrs_enY = []; qrs_enZ = [];
if isempty(Fid_pts_X)
    qrsonX = [];
    qrsoffX = [];
    PACfeatures.maxRRi = NaN;
else
    %dodalem tu wyznaczanie RRow wg algorytmu Kasi Muter - on uwzglednia NaNy
    %pomiedzy Rami - matsol 09042020
    
    %poki co uzywam tylko R z odprowadzenia X - bedzie trzeba pomyslec nad
    %wyznaczeniem Rów globalnie
    for i = 1:length(Fid_pts_X.R)-1
        if (isnan(Fid_pts_X.R(i))==true || isnan(Fid_pts_X.R(i+1))==true)
            continue;
        else
            RR_temp = Fid_pts_X.R(i+1) - Fid_pts_X.R(i);
        end
        RR_interval = [RR_interval, RR_temp];
    end
    PACfeatures.maxRRi = max(RR_interval)/Fs;
    
    qrsonX = Fid_pts_X.QRSon(~(isnan(Fid_pts_X.QRSon) | isnan(Fid_pts_X.QRSoff)));
    qrsoffX = Fid_pts_X.QRSoff(~(isnan(Fid_pts_X.QRSon) | isnan(Fid_pts_X.QRSoff)));
    for i = 1:length(qrsonX)
        qrs_fragX = sigX(qrsonX(i):qrsoffX(i));
        qrs_enX = [qrs_enX; trapz((qrs_fragX).*(qrs_fragX))];
    end
end

if isempty(Fid_pts_Y)
    qrsonY = [];
    qrsoffY = [];
else
    qrsonY = Fid_pts_Y.QRSon(~(isnan(Fid_pts_Y.QRSon) | isnan(Fid_pts_Y.QRSoff)));
    qrsoffY = Fid_pts_Y.QRSoff(~(isnan(Fid_pts_Y.QRSon) | isnan(Fid_pts_Y.QRSoff)));
    for i = 1:length(qrsonY)
        qrs_fragY = sigY(qrsonY(i):qrsoffY(i));
        qrs_enY = [qrs_enY; trapz((qrs_fragY).*(qrs_fragY))];
    end
end

if isempty(Fid_pts_Z)
    qrsonZ = [];
    qrsoffZ = [];
else
    qrsonZ = Fid_pts_Z.QRSon(~(isnan(Fid_pts_Z.QRSon) | isnan(Fid_pts_Z.QRSoff)));
    qrsoffZ = Fid_pts_Z.QRSoff(~(isnan(Fid_pts_Z.QRSon) | isnan(Fid_pts_Z.QRSoff)));
    for i = 1:length(qrsonZ)
        qrs_fragZ = sigZ(qrsonZ(i):qrsoffZ(i));
        qrs_enZ = [qrs_enZ; trapz((qrs_fragZ).*(qrs_fragZ))];
    end
end

PVCfeatures.qrs_en_std = max( [std(qrs_enX), std(qrs_enY), std(qrs_enZ)])/Fs;
PVCfeatures.qrs_en_max_min = max([max(qrs_enX) - min(qrs_enX), max(qrs_enY) ...
    - min(qrs_enY), max(qrs_enZ) - min(qrs_enZ)])/Fs;

qrs_lenX = qrsoffX - qrsonX;
qrs_lenY = qrsoffY - qrsonY;
qrs_lenZ = qrsoffZ - qrsonZ;
PVCfeatures.qrs_len_std = max([std(qrs_lenX), std(qrs_lenY), std(qrs_lenZ)])/Fs;
PVCfeatures.qrs_len_max_min = max([max(qrs_lenX) - min(qrs_lenX), max(qrs_lenY) ...
    - min(qrs_lenY), max(qrs_lenZ) - min(qrs_lenZ)])/Fs;

if isempty(Fid_pts_2chann)
    PVCfeatures.t_ampl_anyneg = NaN;
    PVCfeatures.qrsoff_anyneg = NaN;
else
    t_ampl = sig2(Fid_pts_2chann.T(~isnan(Fid_pts_2chann.T)));
    PVCfeatures.t_ampl_anyneg = any(t_ampl<0);
    
    sig2_qrsoff = sig2(Fid_pts_2chann.QRSoff(~isnan(Fid_pts_2chann.QRSoff)));
    PVCfeatures.qrsoff_anyneg = any(diff(sig2_qrsoff)<0);
end
% RR_i+1>RR_i & (RR_i-1 >= 0.90*(RR_i+1+RR_i)/2 ||  RR_i-1<= 1.10*(RR_i+1+RR_i)/2)
RR = RR_interval;
%zmienilem tutaj kod na nowy - matsol 09042020
if(length(RR)>2)
    PVCfeatures.anypause = any(RR(3:end)>RR(2:end-1) & (RR(1:end-2) >= (0.90*(RR(3:end)+RR(2:end-1))/2) ...
        |  RR(1:end-2)<= (1.10*(RR(3:end)+RR(2:end-1))/2)));
else
    PVCfeatures.anypause = 0;
end
%
% RRi = Fid_pts_XYZ.R';
% RRi(1) = []; RRi(end) = [];
% RRiplus1 = Fid_pts_XYZ.R';
% RRiplus1(1:2) = [];
% RRiminus1 = Fid_pts_XYZ.R';
% RRiminus1(end-1:end) = [];
% RRi = diff(RRi);
% RRiplus1 = diff(RRiplus1);
% RRiminus1 = diff(RRiminus1);
% PVCfeatures.anypause = any(RRiplus1>RRi & (RRiminus1 >= ...
%     (0.90*(RRiplus1+RRi)/2) |  RRiminus1<= (1.10*(RRiplus1+RRi)/2)));
end

function Fid_pts_sig = fiducial_points(ECG, header_data)
% input:
% ECG - ecg signal (one lead signal) for example ECG12filt(1,:)
% header_data
% output:
% fiducial points for one lead (that lead which was given as ecg
% signal)

addpath(genpath('Tools/'))
load('HRVparams_12ECG','HRVparams')
[recording,Total_time,~,Fs,~,~,~]=extract_data_from_header(header_data);

HRVparams.Fs=Fs;
HRVparams.PeakDetect.windows = floor(Total_time-1);
HRVparams.windowlength = floor(Total_time);

sig = [ECG;ECG;ECG]';
VecMag = vecnorm(sig');

[~, ~, jqrs_ann, ~ , ~] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
%sqi = [tSQI', SQIvalue'];
ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
wavedet_config.setup.wavedet.QRS_detection_only = 0;
[Fid_pts_sig,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
end

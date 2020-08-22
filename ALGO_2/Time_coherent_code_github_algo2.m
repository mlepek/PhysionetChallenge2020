% Input:
%      1. sampling frequency in HZ
%      2. matlab file containing the following:
%        - XYZ_O(:,3): 10 seconds X, Y, and Z leads
%        - Rx: index of R peak on X lead
%        - Ry: index of R peak on Y lead
%        - Rz: index of R peak on Z lead

% Output:
%      2. Mat file containing the calculated median beat



function [XYZ_M_T,Fid_pts_Median] = Time_coherent_code_github_algo2(XYZLeads,Fid_pts,fs)


%% ===================== Initialize variables ========================
Rpeaks=[];
num_leads = size(XYZLeads, 2);
XYZ_M_T=zeros(120*(fs/100), num_leads);


%% ===================== load variables from .mat file ========================
XYZ_O=XYZLeads;



Rx = Fid_pts.qrs;
Rx(Rx<200)=[];
Rx(Rx>(length(XYZLeads)-200))=[];
Rpeaks = Rx';

if length(Rx) < 2
    Fid_pts_Median = [];
    return
end

R_no=length(Rx);

%% Total number of samples
Total_samples=length(XYZ_O(:,1));

%% Length of each beat is set to 1200 ms.
Beat_length=fs*1.2;
BLengthD2=floor(Beat_length/2);

% =========================== Variable Calculation ============-================
if ~isempty(Rpeaks)
c_int_init=1;
while (Rpeaks(c_int_init,1)<(BLengthD2))
        c_int_init=c_int_init+1;
end

%% Exclude the last if we do not have a complete last beat
%Poprawka, bo czêsto Rpeaks ma tylko jeden row czego ta funkcja nie
%obs³ugiwa³a
if size(Rpeaks,1)>1 
    if abs(Rpeaks(end-1,1)-Total_samples)<(BLengthD2) 
        c_fin=2;
    else
        c_fin=1;
    end
else 
    if abs(Rpeaks(1)-Total_samples)<(BLengthD2)
        c_fin=2;
    else
        c_fin=1;
    end
end
        

for lead_idx = 1:num_leads
    c_int = c_int_init;
    Beats_T = [];
    %% Calculate the absolute maximum first derivative adjacent to each R peak on the each lead.
    I_dv=zeros(R_no,1);

    for ii=c_int:R_no-c_fin
       
        [~,I_dv_temp]=max(abs(diff(XYZ_O(Rx(ii)-199:Rx(ii)+200,lead_idx))));
        I_dv(ii)=Rx(ii)-(199-I_dv_temp);
       
    end


    %% Separate the each beat in each lead. Each beat is centered on Maximum dv/dt detected.
    %To samo trzeba by³o uwzglêdniæ w tym miejscu
    if R_no>1
        while (I_dv(c_int)<BLengthD2)
                c_int=c_int+1;
        end
    end


    bct=1;
    for ii=c_int:R_no-c_fin
        if I_dv(ii)+(BLengthD2) > size(XYZ_O,1)
            break;
        end

        Beats_T(:,bct)=XYZ_O(I_dv(ii)-((BLengthD2)-1):I_dv(ii)+(BLengthD2),lead_idx);

        bct=bct+1;
    end
    
    %% Create the median beat of each lead
    XYZ_M_T(:,lead_idx) = median(Beats_T');
end

end
%% Calculate Median Fid points

QRSon_diff = nanmedian(Fid_pts.qrs - Fid_pts.QRSon);
QRSoff_diff = nanmedian(Fid_pts.QRSoff - Fid_pts.qrs);
Tpeak_diff = nanmedian(Fid_pts.T - Fid_pts.qrs);
Toff_diff = nanmedian(Fid_pts.Toff - Fid_pts.qrs);

Fid_pts_Median.QRS = floor(Beat_length/2);

Fid_pts_Median.QRSon = floor(Fid_pts_Median.QRS - QRSon_diff);
Fid_pts_Median.QRSoff = floor(Fid_pts_Median.QRS + QRSoff_diff);
Fid_pts_Median.Tpeak = floor(Fid_pts_Median.QRS + Tpeak_diff);
Fid_pts_Median.Toff = floor(Fid_pts_Median.QRS + Toff_diff);


end


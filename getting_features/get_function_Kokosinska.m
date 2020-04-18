function [TFRPrimaV1V2V3,RSmax,QRSwidth,SyRyMedianV1,QRpercmed,QRStimemin] = get_function_Kokosinska(ECG12filt,XYZLeads,header_data)
%RBBB:
%TFRPrimaV1V2V3 - czy R' wystêpuje w którymkolwiek z odprowadzeñ V1, V2 lub V3
%RSmax-Maksumalna wartoœæ amplitudy za³amka S do R z odprowadzeñ I, aVL, V5 i V6
%LBBB:
%QRSwidth - d³ugoœæ zespo³u QRS  (z uœrednionego beatu)
%SyRyMedianV1 -stosunek amplitud za³amków S i R w odprowadzeniu V1
%QRpercmed- procent liczby za³amków Q w stosunku do liczby za³amków R - mediana z wartoœci obliczonych w odprowadzeniach I, V5 i V6 
%QRStimemin - czas trwania od pocz¹ku QRS do peaku R w odprowadzeniach V5 & V6 - najmniejsza wartoœæ uzyskan¹ z uœrednionych beatów z odprowadzeñ V5 i V6

addpath(genpath('Tools/'))
load('HRVparams_12ECG','HRVparams')

[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

HRVparams.Fs=Fs;
HRVparams.PeakDetect.windows = floor(Total_time-1);
HRVparams.windowlength = floor(Total_time);
VecMag = vecnorm(XYZLeads');



%---------------------------RBBB-------------------------------------------

%TFRPrimaV1V2V3------------------------------------------------------------

%Œredni beat w odprowadzeniach odpowiednio V1, V2 i V3
leadsV1V2V3=[7 8 9];%7 8 9
for j=1:length(leadsV1V2V3)
    i = leadsV1V2V3(j); 
%     disp(i)
    V1V2V3 = [ECG12filt(i,:);ECG12filt(i,:);ECG12filt(i,:)]';
    VecMag = vecnorm(V1V2V3');
    try
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
         sqi = [tSQI', SQIvalue'];
         ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        [Fid_pts_V1V2V3,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
    %   [XYZ_Median_V1V2V3,Fid_pts_Median_V1V2V3] = Time_coherent_code_github(V1V2V3,Fid_pts_V1V2V3,Fs);

        %Rprima 
        Rprima = [Fid_pts_V1V2V3(:).Rprima];
        %Finding if ANY NaN values occur in a matrix.
        ResultRprima = any(isnan(Rprima(:)));

        %TFRPrima - czy R' wystêpuje w odprowadzeniu V1/V2/V3(1-tak,0-nie)
        TFRPrima(j)=uint8(~ResultRprima);
    catch
        TFRPrima(j)=0;
    end
   

end

%TFRPrimaV1V2V3 - czy R' wystêpuje w którymkolwiek z odprowadzeñ V1, V2 lub V3
if any(TFRPrima>0)
    TFRPrimaV1V2V3=1;
else
    TFRPrimaV1V2V3=0;
end


%RSmax---------------------------------------------------------------------

%Odprowadzenia I & aVL & V5 & V6 i obliczenia œredniego beatu 
leadsIaVLV5V6=[1 11 12]; % nie uwzglêdni³am aVL (5), bo coœ tam jest nie tak :(
for j=1:length(leadsIaVLV5V6)
    i = leadsIaVLV5V6(j);
    IaVLV5V6 = [ECG12filt(i,:);ECG12filt(i,:);ECG12filt(i,:)]';
    VecMag = vecnorm(IaVLV5V6');
    try
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
        sqi = [tSQI', SQIvalue'];
        ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        [Fid_pts_IaVLV5V6,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
        [XYZ_Median_IaVLV5V6,Fid_pts_Median_IaVLV5V6] = Time_coherent_code_github(IaVLV5V6,Fid_pts_IaVLV5V6,Fs);
        QRStart=[Fid_pts_Median_IaVLV5V6(:).QRSon];
        QRSend=[Fid_pts_Median_IaVLV5V6(:).QRSoff];
        %Za³amek R 
        if all(isnan(XYZ_Median_IaVLV5V6))==1
            SyRyMedian(j)=0;
        else if max(XYZ_Median_IaVLV5V6)-min(XYZ_Median_IaVLV5V6)<0.05
            SyRyMedian(j)=0;
        else 
            values=XYZ_Median_IaVLV5V6(QRStart:QRSend);
            %Amplituda R
            Ry=max(values);
            %Amplituda S
            Sy=min(values);
            SyRyMedian(j)=Sy/Ry;
            end
        end
    catch
        SyRyMedian(j)=0;
    end
end
if SyRyMedian ~=0  
    %RSmax-Maksumalna wartoœæ amplitudy za³amka S do R z odprowadzeñ I, aVL, V5 i V6
    RSmax=max(SyRyMedian);
else
    RSmax=0;
end



%---------------------------LBBB-------------------------------------------

%QRSwidth------------------------------------------------------------------
VecMag = vecnorm(XYZLeads');
try
    % Convert ECG waveform in rr intervals
    [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
    sqi = [tSQI', SQIvalue'];

    % Find fiducial points using ECGKit
    ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
    wavedet_config.setup.wavedet.QRS_detection_only = 0;
    [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
    [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);
    QRSwidth = [Fid_pts_Median(:).QRSoff]-[Fid_pts_Median(:).QRSon];
catch
    QRSwidth = 0;
end

%SyRyMedianV1--------------------------------------------------------------

%Odprowadzenie V1 i obliczenia œredniego beatu 
V1 = [ECG12filt(7,:);ECG12filt(7,:);ECG12filt(7,:)]';
VecMag = vecnorm(V1');
try
    [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
    %sqi = [tSQI', SQIvalue'];
    ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
    wavedet_config.setup.wavedet.QRS_detection_only = 0;
    [Fid_pts_V1,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
    [XYZ_Median_V1,Fid_pts_Median_V1] = Time_coherent_code_github(V1,Fid_pts_V1,Fs);
    QRStart=[Fid_pts_Median_V1(:).QRSon];
    QRSend=[Fid_pts_Median_V1(:).QRSoff];
        if all(isnan(XYZ_Median_V1))==1
            SyRyMedianV1=0;
        else if max(XYZ_Median_V1)-min(XYZ_Median_V1)<0.05
            SyRyMedianV1=0;
        else
            %Za³amek R 
            values=XYZ_Median_V1(QRStart:QRSend);
            %Amplituda R
            Ry=max(values);
            %Amplituda S
            Sy=min(values);
            %SyRyMedianV1 -stosunek amplitud za³amków S i R w odprowadzeniu V1
            SyRyMedianV1=Sy/Ry;
            end
        end
catch
    SyRyMedianV1=0;
end

%QRperc--------------------------------------------------------------------

%Odprowadzenia I & V5 & V6 i obliczenia œredniego beatu 
leadsIV5V6=[1 11 12];
for j=1:length(leadsIV5V6)
    i = leadsIV5V6(j);
    IV5V6 = [ECG12filt(i,:);ECG12filt(i,:);ECG12filt(i,:)]';
    VecMag = vecnorm(IV5V6');
    try
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
       % sqi = [tSQI', SQIvalue'];
        ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        [Fid_pts_IV5V6,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
    %     [XYZ_Median_IV5V6,Fid_pts_Median_IV5V6] = Time_coherent_code_github(IV5V6,Fid_pts_IV5V6,Fs);
        %numR - liczba za³amków R w odprowadzeniach I & V5 & V6
        numR=numel([Fid_pts_IV5V6(:).R]);
        %number of non nan observations
        numQ=sum(~isnan(([Fid_pts_IV5V6(:).Q])));
        %QRperc- procent liczby za³amków Q w stosunku do liczby za³amków R
        QRperc(j)=100*numQ/numR;
    catch
       QRperc(j)=0;
    end
end

%QRpercmed - mediana z wartoœci obliczonych w odprowadzzeniach I, V5 i V6
QRpercmed=nanmedian(QRperc);
%QRStime-------------------------------------------------------------------


%Odprowadzenia V5 & V6 i obliczenia œredniego beatu 
leadsV5V6=[11 12];
for j=1:length(leadsV5V6)
    i = leadsV5V6(j);
    V5V6 = [ECG12filt(i,:);ECG12filt(i,:);ECG12filt(i,:)]';
    VecMag = vecnorm(V5V6');
    try
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
        %sqi = [tSQI', SQIvalue'];
        ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        [Fid_pts_V5V6,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);
        [XYZ_Median_V5V6,Fid_pts_Median_V5V6] = Time_coherent_code_github(V5V6,Fid_pts_V5V6,Fs);
        %QRStime - czas trwania od pocz¹ku QRS do peaku R w odprowadzeniach V5 & V6
        QRStart=[Fid_pts_Median_V5V6(:).QRSon];
        QRSend=[Fid_pts_Median_V5V6(:).QRSoff];
        if all(isnan(XYZ_Median_V5V6))==1
            QRStime(j)=0;
        else if max(XYZ_Median_V5V6)-min(XYZ_Median_V5V6)<0.05
            QRStime(j)=0;  
        else
            %Za³amek R 
            values=XYZ_Median_V5V6(QRStart:QRSend);
            Ry=max(values);
            Rx=find(XYZ_Median_V5V6==Ry);
            QRStime(j)=Rx(1,:)-QRStart;    
            end
        end
    catch
       QRStime(j)=0;
    end
QRStimemin=min(QRStime);
end


function   coherentMatrixToExport = getAveragedBeat(fs_fixed, data,ann, max_length)

   addpath('./ALGO_2/');
   addpath(genpath('./Tools/'));

    coherentMatrixToExport = zeros(12, 1.2*fs_fixed); 
    wavedet_config.setup.wavedet.QRS_detection_only = 0;
    ECG_header.nsig = 1; ECG_header.freq = fs_fixed; 
%    
matrixToExport_withzeros = zeros(12,max_length);
    if(size(data,2)>size(matrixToExport_withzeros,2)) %ucinanie dluzszych sygnalow
        data = data(:,1:size(matrixToExport_withzeros,2));
    end
%     
    ECG_header.nsamp = length(data);
%     ann = hierarchiczny_12(data',fs_fixed);
    try
        [Fid_pts,~,~] = wavedet_3D_ECGKit(data(1,:)', ann', ECG_header, wavedet_config);
        [coherent,~] = Time_coherent_code_github_algo2(data',Fid_pts,fs_fixed);
        coherentMatrixToExport(:,1:length(coherent)) = coherent';
    catch
    end
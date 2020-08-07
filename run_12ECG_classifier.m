function [score, label,classes] = run_12ECG_classifier(data,header_data, loaded_model)

        %tline = fgetl(header_data);
        %tmp_str = strsplit(tline,' ');
        tmp_str = strsplit(header_data{1}, ' ');
%         recording_label = tmp_str{1};   
        fs = str2num(tmp_str{3});


    dlnet1 = loaded_model.model;
	%model = loaded_model.model;
	classes = loaded_model.classes;
    num_classes = length(classes);
    %label = zeros([1,num_classes]);
    %score = ones([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    %features = get_12ECG_features(data,header_data);

    %preprocessing
    
    %PARAMETRY SYGNALOW
    fs_fixed = 100; %docelowe sample frequency
    max_length = round(10*fs_fixed); %ustawienie maksymalnej dlugosci sygnalu: 10 s%DOWNSAMPLING
    data= resample(data',fs_fixed,fs); %resampling do 100Hz
    data = data';   
    
    %DATA FILTERING
    ECG12filt = [];
    % median filter to remove bw
    for ii=1:12
        ECG12filt(ii,:) = movmean(data(ii,:)', 5);
    end
    ECG12filt2 = ECG12filt;
    
    for iii=1:size(ECG12filt,1)
        ECG12filt2(iii,:) = ECG12filt(iii,:) - movmedian(ECG12filt(iii,:), 50);
    end
    data = ECG12filt2;   
    
    %PRZYCINANIE DLUZSZYCH SYGNALOW
    matrixToExport_withzeros = zeros(12,max_length);
    if(size(data,2)>size(matrixToExport_withzeros,2)) %usinanie dluzszych sygnalow
        matrixToExport_withzeros(1:12,1:size(matrixToExport_withzeros,2)) = data(:,1:size(matrixToExport_withzeros,2));
    else
        matrixToExport_withzeros(1:12,1:size(data,2)) = data;
    end
    %             all_signals_matrix = [all_signals_matrix reshape(matrixToExport_withzeros',[],1)];
    signal_to_predict = reshape(matrixToExport_withzeros',[],1); 
    
    X_Train = signal_to_predict;
    
    %przegrupowanie macierzy do 3d
    X_Train = reshape(X_Train,size(X_Train,1)/12,12,size(X_Train,2));%przegrupowanie macierzy do 4d
    X_Train_Input = zeros(size(X_Train,1),size(X_Train,2),1,size(X_Train,3),'int16');
    X_Train_Input(:,:,1,1) = X_Train(:,:,1);
    
    XTrain1 = X_Train_Input;
    
    dlXTrain1 = dlarray(single(XTrain1),'SSCB');
    %if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
    %    dlXTrain1 = gpuArray(dlXTrain1);
    %end
    
    labelThreshold = 0.5;
    
    dlYPred1 = forward(dlnet1,dlXTrain1);
    dlYPred_concat_train = sigmoid(dlYPred1);
    YPred_Train_toPerformance = extractdata(dlYPred_concat_train> labelThreshold);
    
    score = extractdata(dlYPred_concat_train);
    label =  YPred_Train_toPerformance;
    
    classes = cellstr(num2str(classes)); % to potrzebne bo driver wymaga tablicy cell
    
    %score = mnrval(model,features);		
    %[~,idx] = max (score);
    %label(idx)=1;
end




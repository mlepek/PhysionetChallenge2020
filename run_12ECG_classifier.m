function [score, label,classes] = run_12ECG_classifier(data, header_data, loaded_model)

    %persistent mp;
    %persistent wyniki_scores;
    %persistent wyniki_labels;

        %tline = fgetl(header_data);
        %tmp_str = strsplit(tline,' ');
        tmp_str = strsplit(header_data{1}, ' ');
%         recording_label = tmp_str{1};   
        fs = str2num(tmp_str{3});


    dlnet1 = loaded_model.model1;
    dlnet2 = loaded_model.model2;
    dlnet3 = loaded_model.model3;
    dlnet = loaded_model.model;
    
	%model = loaded_model.model;
	classes = loaded_model.classes;
    num_classes = length(classes);
    %label = zeros([1,num_classes]);
    %score = ones([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    %features = get_12ECG_features(data,header_data);

    %% PRZYGOTOWANIE SYGNALOW ECG 
    
    %PARAMETRY SYGNALOW
     fs_fixed = 100; %docelowe sample frequency
    max_length = round(30*fs_fixed); %ustawienie maksymalnej dlugosci sygnalu: 10 s %DOWNSAMPLING
%     data= resample(data',fs_fixed,fs); %resampling do 100Hz
%     data = data';   
        
    signal_to_predict = preprocessing_before_aggregation(data, fs, max_length);

%     %DATA FILTERING
%     ECG12filt = [];
%     % median filter to remove bw
%     for ii=1:12
%         ECG12filt(ii,:) = movmean(data(ii,:)', 5);
%     end
%     ECG12filt2 = ECG12filt;
%     
%     for iii=1:size(ECG12filt,1)
%         ECG12filt2(iii,:) = ECG12filt(iii,:) - movmedian(ECG12filt(iii,:), 50);
%     end
%     data = ECG12filt2;   
%     
%     data_new = [];
%     
%     %PRZYCINANIE DLUZSZYCH SYGNALOW
%     matrixToExport_withzeros = zeros(12,max_length);
% %     if(size(data,2)>size(matrixToExport_withzeros,2)) %usinanie dluzszych sygnalow
% %         matrixToExport_withzeros(1:12,1:size(matrixToExport_withzeros,2)) = data(:,1:size(matrixToExport_withzeros,2));
% %     else
% %         matrixToExport_withzeros(1:12,1:size(data,2)) = data;
% %     end
% %     %             all_signals_matrix = [all_signals_matrix reshape(matrixToExport_withzeros',[],1)];
% 
%      if(size(data,2)>size(matrixToExport_withzeros,2)) %usinanie dluzszych sygnalow
%         matrixToExport_withzeros(1:12,1:size(matrixToExport_withzeros,2)) = data(:,1:size(matrixToExport_withzeros,2));
%      else
%         %kopiowanie za krotkich sygnalow:
%         if(size(data,2)<size(matrixToExport_withzeros,2))
%             howManyTimes = floor(size(matrixToExport_withzeros,2)/size(data,2)); %jaka jest wielokrotnosc dlugosci sygnalu (zaokraglona w dol)
%             howManyToCopy_residuum = size(matrixToExport_withzeros,2)-size(data,2)*howManyTimes; %reszta brakujacych probek
%             for j=1:howManyTimes
%                 data_new = [data_new data];
%             end
%            
%             %jesli zostaly jakies resztki (po dodaniu wielokrotnosci)
%             if(howManyToCopy_residuum>0)
%                 data_new = [data_new data(:,howManyToCopy_residuum)];
%             end
%             data = data_new;   
%         end
%        
%         matrixToExport_withzeros(1:12,1:size(data,2)) = data;
%     end
% 
%     signal_to_predict = reshape(matrixToExport_withzeros',[],1); 
%     
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
    
    
    %% PRZYGOTOWANIE SYGNALU RR
    
    [RRsignal,ann] = getRR(max_length,fs_fixed, data);
    dlRRsignal = dlarray(single(RRsignal),'SSCB');
    
    %% PRZYGOTOWANIE USREDNIONYCH BEATOW
    averagedBeat = getAveragedBeat(fs_fixed, data,ann, max_length);
    dlaveragedBeat = dlarray(single(averagedBeat),'SSCB');
    
    %% KLASYFIKACJA
    numNets = 3;
    numHiddenDimension = 20;
    labelThreshold = 0.35;
    
        dlY_Train_Pred1 = predict(dlnet1,dlXTrain1);
        dlY_Train_Pred2 = predict(dlnet2,dlaveragedBeat);
        dlY_Train_Pred3 = predict(dlnet3,dlRRsignal);
        dlX_concat=[dlY_Train_Pred1;dlY_Train_Pred2;dlY_Train_Pred3];
        dlX_concat=reshape(dlX_concat,[1 numNets*numHiddenDimension, 1, size(dlX_concat, 2)]);
        dlXTrain=dlarray(single(dlX_concat),'SSCB');
        dlYPred=gather(predict(dlnet,dlXTrain));
    
    %dlYPred1 = forward(dlnet1,dlXTrain1);
%     dlYPred1 = predict(dlnet1,dlXTrain1);
    dlYPred_concat_train = sigmoid(dlYPred);
    YPred_Train_toPerformance = extractdata(dlYPred_concat_train> labelThreshold);
    
    score = extractdata(dlYPred_concat_train);
    label =  YPred_Train_toPerformance;
    
    % upewnienie sie ze scores i label to wektory wierszowe
    score = score(:)';
    label = label(:)';
    
    %classes = cellstr(num2str(classes)); % to potrzebne bo driver wymaga tablicy cell
   
    %score = mnrval(model,features);		
    %[~,idx] = max (score);
    %label(idx)=1;
    
    % here we do some things to follow the format with 111 classes
    
    all_111_classes = [10370003,111288001,11157007,111975006,13640000,164861001,164865005,164867002,164873001,164884008,164889003,164890007,164895002,164896001,164909002,164917005,164921003,164930006,164931005,164934002,164937009,164947007,164951009,17338001,195042002,195060002,195080001,195101003,195126007,204384007,233917008,251120003,251139008,251146004,251164006,251168009,251170000,251173003,251180001,251182009,251200008,251259000,251266004,251268003,253339007,253352002,266249003,266257000,270492004,27885002,282825002,284470004,29320008,314208002,368009,370365005,39732003,413444003,413844008,425419005,425623009,425856008,426177001,426434006,426627000,426648003,426664006,426749004,426761007,426783006,426995002,427084000,427172004,427393009,428417006,428750005,429622005,445118002,445211001,446358003,446813000,47665007,49260003,49578007,53741008,54016002,54329005,55930002,57054005,59118001,59931005,60423000,63593006,6374002,65778007,67198005,67741000119109,698247007,698252002,704997005,713422000,713426002,713427006,74390002,74615001,75532003,77867006,81898007,82226007,84114007,89792004];
    all_111_classes = all_111_classes';
    %all_111_classes = cellstr(num2str(all_111_classes));
    
    all_111_scores = zeros(1,111);
    all_111_labels = zeros(1,111);
    
    % teraz nadaje niektorym polom w all_111 wartosci (oceniane klasy)
    
    for i = 1:length(classes)
    
        position_of_this_class = find(all_111_classes == classes(i)); 

        all_111_scores(position_of_this_class) = score(i);
        all_111_labels(position_of_this_class) = label(i);
      
    end
    
    score = all_111_scores;
    label = all_111_labels;
    classes = all_111_classes;
    classes = cellstr(num2str(classes));
    classes = strtrim(classes);
    
    
    %mp = [mp; label];
    %wyniki_scores = [wyniki_scores; score];
    %wyniki_labels = [wyniki_labels; label];
    
end




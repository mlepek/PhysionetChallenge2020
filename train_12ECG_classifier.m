function  model = train_12ECG_classifier(input_directory,output_directory)

disp('Loading data...')

% Find files.
input_files = {};
for f = dir(input_directory)'
    if exist(fullfile(input_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'mat')
        input_files{end + 1} = f.name;
    end
end

% read number of unique classes
 %classes_to_save = get_classes(input_directory,input_files); 
%classes_file = importfile('./dx_mapping_scored.csv'); % to sie
%wywrocilo podczas processingu entry
%classes = classes_file.SNOMEDCTCode; %kody klas
classes = [270492004;
           164889003;
           164890007;
           426627000;
           713427006;
           713426002;
           445118002;
            39732003;
           164909002;
           251146004;
           698252002;
            10370003;
           284470004;
           427172004;
           164947007;
           111975006;
           164917005;
            47665007;
            59118001;
           427393009;
           426177001;
           426783006;
           427084000;
            63593006;
           164934002;
            59931005;
            17338001];
% classes = {};
% for i=1:27
%     classes{i} = classes_temp(i);
% end

num_classes = length(classes);
num_files = length(input_files);
Total_data=cell(1,num_files);
Total_header=cell(1,num_files);
Labels_hotmatrix = zeros(num_classes,num_files);
%PARAMETRY SYGNALOW
fs_fixed = 100; %docelowe sample frequency
max_length = round(10*fs_fixed); %ustawienie maksymalnej dlugosci sygnalu: 10 s

% Iterate over files.
all_signals_matrix = zeros(12*max_length,num_files,'int16'); %inicjalizacja macierzy z sygnalami do eksportu
iter=1;

for i = 1:num_files
    disp(['    ', num2str(i), '/', num2str(num_files), '...'])
    
    % Load data.
    file_tmp=strsplit(input_files{i},'.');
    tmp_input_file = fullfile(input_directory, file_tmp{1});
    
    [~,~,single_recording_labels, fs]=get_true_labels([tmp_input_file '.hea'],classes);
    
    [data,hea_data] = load_challenge_data(tmp_input_file);
    if(sum(single_recording_labels==1)>0)
        %Preprocessing of the data before aggregation to one matrix
        all_signals_matrix = preprocessing_before_aggregation(data, fs, all_signals_matrix, iter);
        
        
        Labels_hotmatrix(:,iter) = single_recording_labels';
        
        Total_data{i}=data;
        Total_header{i}=hea_data;
        iter = iter + 1;
    end
    
end

all_signals_matrix = all_signals_matrix(:,1:iter);
Labels_hotmatrix = Labels_hotmatrix(:,1:iter);
disp('Training model..')

% label=zeros(num_files,num_classes);
%
% for i = 1:num_files
%
%     disp(['    ', num2str(i), '/', num2str(num_files), '...']);
%
%     data = Total_data{i};
%     header_data = Total_header{i};
%
%     tmp_features = get_12ECG_features(data,header_data);
%
%     features(i,:)=tmp_features;
%
%     for j = 1 : length(header_data)
%         if startsWith(header_data{j},'#Dx')
%             tmp = strsplit(header_data{j},': ');
%             tmp_c = strsplit(tmp{2},',');
%             for k=1:length(tmp_c)
%                 idx=find(strcmp(classes,tmp_c{k}));
%                 label(i,idx)=1;
%             end
%             break
%         end
%     end
%
%
% end


X_Train = all_signals_matrix;
% Labels_Train = Labels_rand;
% Labels_abbr_Train = Labels_abbr_rand;
Labels_hotmatrix_Train = Labels_hotmatrix;

Y_Train = Labels_hotmatrix_Train;
%przegrupowanie macierzy do 3d
X_Train = reshape(X_Train,size(X_Train,1)/12,12,size(X_Train,2));
%przegrupowanie macierzy do 4d
X_Train_Input = zeros(size(X_Train,1),size(X_Train,2),1,size(X_Train,3),'int16');
for i=1:length(X_Train)
    X_Train_Input(:,:,1,i) = X_Train(:,:,i);
end

XTrain1 = X_Train_Input;
YTrain = Y_Train;

%czyszczenie pamieci
clear X_Train

%% MIEJSCE NA TRAINING CUSTOM LOOP

layers = [
    imageInputLayer([size(X_Train_Input,1) size(X_Train_Input,2) 1],"Name","imageinput","Mean",mean(X_Train_Input,4))
    
    convolution2dLayer([3,3],8,'Padding','same',"Name","conv_1")
    batchNormalizationLayer("Name","batchnorm_1")
    reluLayer("Name","relu_1")
    convolution2dLayer([3,3],8,'Padding','same',"Name","conv_2")
    batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_2")
    maxPooling2dLayer([2,2],"Name","maxpool_1",'Stride',[2,2],'Padding','same')
    
    %dropoutLayer(0.1,"Name","dropout_1")
    
    convolution2dLayer(3,16,'Padding','same',"Name","conv_3")
    batchNormalizationLayer("Name","batchnorm_3")
    reluLayer("Name","relu_3")
    convolution2dLayer(3,16,'Padding','same',"Name","conv_4")
    batchNormalizationLayer("Name","batchnorm_4")
    reluLayer("Name","relu_4")
    maxPooling2dLayer(2,'Stride',2,"Name","maxpool_2",'Padding','same')
    
    %dropoutLayer(0.1,"Name","dropout_2")
    
    convolution2dLayer(3,32,'Padding','same',"Name","conv_5")
    batchNormalizationLayer("Name","batchnorm_5")
    reluLayer("Name","relu_5")
    convolution2dLayer(3,32,'Padding','same',"Name","conv_6")
    batchNormalizationLayer("Name","batchnorm_6")
    reluLayer("Name","relu_6")
    maxPooling2dLayer(2,'Stride',2,"Name","maxpool_3",'Padding','same')
    
    %dropoutLayer(0.1,"Name","dropout_3")
    
    convolution2dLayer(3,64,'Padding','same',"Name","conv_7")
    batchNormalizationLayer("Name","batchnorm_7")
    reluLayer("Name","relu_7")
    convolution2dLayer(3,64,'Padding','same',"Name","conv_8")
    batchNormalizationLayer("Name","batchnorm_8")
    reluLayer("Name","relu_8")
    
    %dropoutLayer(0.1,"Name","dropout_4")
    %{
    convolution2dLayer(3,128,'Padding','same',"Name","conv_9")
    batchNormalizationLayer("Name","batchnorm_9")
    reluLayer("Name","relu_9")
    convolution2dLayer(3,128,'Padding','same',"Name","conv_10")
    batchNormalizationLayer("Name","batchnorm_10")
    reluLayer("Name","relu_10")
    
    dropoutLayer(0.2,"Name","dropout_5")
    
    convolution2dLayer(3,256,'Padding','same',"Name","conv_11")
    batchNormalizationLayer("Name","batchnorm_11")
    reluLayer("Name","relu_11")
    convolution2dLayer(3,256,'Padding','same',"Name","conv_12")
    batchNormalizationLayer("Name","batchnorm_12")
    reluLayer("Name","relu_12")
    
    dropoutLayer(0.2,"Name","dropout_6")
    %}
    fullyConnectedLayer(27,"Name","fully-con")
    ];

lgraph = layerGraph(layers);
dlnet1 = dlnetwork(lgraph);


%% Specify the training options

numEpochs = 1;
miniBatchSize = 1024;
epsilon=0.001;
learnRate = 0.001;
GradDecay=0.9;
sqGradDecay= 0.9;
executionEnvironment = "cpu";
numHiddenDimension = 27;
l2Regularization = 0.06; % 0.06 wypracowana z poprzedniej fazy

labelThreshold = 0.5;

velocity1 = [];
numObservations = size(YTrain,2);
averageSqGrad1=[];
averageGrad1=[];

%plots = "training-progress";

numIterationsPerEpoch = floor(numObservations/miniBatchSize);
validationFrequency = 3 * numIterationsPerEpoch;

%{
if plots == "training-progress"
    figure
    
    % Labeling F-Score.
    subplot(2,1,1)
    lineFScoreTrain = animatedline('Color',[0 0.447 0.741]);
    lineFScoreValidation = animatedline( ...
        'LineStyle','--', ...
        'Marker','o', ...
        'MarkerFaceColor','black');
    ylim([0 1])
    xlabel("Iteration")
    ylabel("Labeling F-Score")
    grid on
    grid(gca,'minor')
    
    % Loss.
    subplot(2,1,2)
    lineLossTrain = animatedline('Color',[0.85 0.325 0.098]);
    lineLossValidation = animatedline( ...
        'LineStyle','--', ...
        'Marker','o', ...
        'MarkerFaceColor','black');
    ylim([0 inf])
    xlabel("Iteration")
    ylabel("Loss")
    grid on
    grid(gca,'minor')
end
%}

iteration = 0;
start = tic;
numClasses = 27;
% Loop over epochs.
for epoch = 1:numEpochs
    %disp(epoch)
    
    % Shuffle data.
    idx = randperm(size(YTrain,2));
    XTrain1 = XTrain1(:,:,:,idx);
    %     XTrain2 = XTrain2(:,:,:,idx);
    YTrain = YTrain(:,idx);
    
    % Loop over mini-batches.
    for i = 1:numIterationsPerEpoch
        iteration = iteration + 1;
        %disp(iteration)
        
        % Read mini-batch of data and convert the labels to dummy
        % variables.
        idx = (i-1)*miniBatchSize+1:i*miniBatchSize;
        X1 = XTrain1(:,:,:,idx);
        
        %X1_Test = XTest1;
        %         X2 = XTrain2(:,:,:,idx);
        % convert the label into one-hot vector to calculate the loss
        %         Y = zeros(numClasses, miniBatchSize, 'single');
        Y = YTrain(:,idx);
        %Y_Test = YTest; --- w jakim to jest celu?
        %         for c = 1:numClasses
        %             Y(c,YTrain(idx)==classes(c)) = 1;
        %             Y = YTrain(:,idx);
        %         end
        
        % Convert mini-batch of data to dlarray.
        dlX1 = dlarray(single(X1),'SSCB');
        %dlX1_Test = dlarray(single(X1_Test),'SSCB'); - w jakim to jest
        %celu tutaj? Przenioslem do czesc walidacyjnej
        %         dlX2 = dlarray(single(X2),'SSCB');
        
        % If training on a GPU, then convert data to gpuArray.
        if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
            dlX1 = gpuArray(dlX1);
            %dlX1_Test = gpuArray(dlX1_Test); -- w jakim to celu?
            %             dlX2 = gpuArray(dlX2);
        end
        %the traning loss and the gradients after the backpropagation were
        %calculated using the helper function modelGradients_demo
        [gradients1,loss,dlYPred_to_draw] = dlfeval(@modelGradients_demo,dlnet1,dlX1,dlarray(Y));
        
        % dodaje regularyzacje L2
        idx = dlnet1.Learnables.Parameter == "Weights";
        gradients1(idx,:) = dlupdate(@(g,w) g + l2Regularization*w, gradients1(idx,:), dlnet1.Learnables(idx,:));
        
        % Update the network parameters using the RMSProp optimizer.
        % Update the parameters in dlnet1 to 3 sequentially
        %         [dlnet3,averageGrad3, averageSqGrad3] = adamupdate(dlnet3, gradients3,averageGrad3, averageSqGrad3,iteration,learnRate,GradDecay,sqGradDecay,epsilon);
        %         [dlnet2,averageGrad2, averageSqGrad2] = adamupdate(dlnet2, gradients2,averageGrad2, averageSqGrad2,iteration,learnRate,GradDecay,sqGradDecay,epsilon);
        [dlnet1,averageGrad1, averageSqGrad1] = adamupdate(dlnet1, gradients1,averageGrad1, averageSqGrad1,iteration,learnRate,GradDecay,sqGradDecay,epsilon);
        % Display the training progress.
        D = duration(0,0,toc(start),'Format','hh:mm:ss');
        %         addpoints(lineLossTrain,iteration,double(gather(extractdata(loss))))
        %         title("Epoch: " + epoch + ", Elapsed: " + string(D))
        %         drawnow
        
        % Display the training progress.
        %{
        if plots == "training-progress"
            subplot(2,1,1)
            D = duration(0,0,toc(start),'Format','hh:mm:ss');
            title("Epoch: " + epoch + ", Elapsed: " + string(D))
            
            % Loss.
            addpoints(lineLossTrain,iteration,double(gather(extractdata(loss))))
            
            % Labeling F-score.
            %             YPred = extractdata(dlYPred_to_draw) > labelThreshold;
            YPred_Train = extractdata(dlYPred_to_draw) == max(dlYPred_to_draw);
            score = labelingFScore(YPred_Train,Y);
            addpoints(lineFScoreTrain,iteration,extractdata(score))
            
            drawnow
            
            %% Display validation metrics.
            %             if iteration == 1 || mod(iteration,validationFrequency) == 0
            %
            %                 % usuwamy dropout, zeby nie dzialal podczas funkcji
            %                 % forward
            %                 % Czy obiekt layerGraph zachowuje gdzies niejawnie
            %                 % informacje o wagach? Wyglada na to, ze tak!
            %                 %
            %                 lgraph = layerGraph(dlnet1);
            %
            %                 larray = [ dropoutLayer(0.0, 'Name','dropout_1') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_1',larray);
            %                 larray = [ dropoutLayer(0.0, 'Name','dropout_2') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_2',larray);
            %                 larray = [ dropoutLayer(0.0, 'Name','dropout_3') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_3',larray);
            %                 larray = [ dropoutLayer(0.0, 'Name','dropout_4') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_4',larray);
            %
            %                 dlnet1 = dlnetwork(lgraph);
            %                 %}
            %
            % %                 [gradients1,lossValidation,dlYPred_to_draw_test] = dlfeval(@modelGradients_demo,dlnet1,dlX1_Test,dlarray(Y_Test));
            %
            %                 %X1_Test = XTest1;
            %                 %Y_Test = YTest;
            %                 dlX1_Test = dlarray(single(XTest1),'SSCB');
            %
            %                 dlX_concat_test=dlarray(single(dlX1_Test),'SSCB');
            %                         %dlX_concat_test = gpuArray(dlX_concat_test);
            %                 dlY_concat_test = forward(dlnet1,dlX_concat_test);
            %                 %dlY_concat_test = predict(dlnet1,dlX_concat_test);
            %                 dlYPred_concat_test = sigmoid(dlY_concat_test);
            %                 lossValidation = crossentropy(dlYPred_concat_test,Y_Test);
            %                 % Loss.
            % %                 lossValidation = crossentropy(dlYPredValidation,TValidation, ...
            % %                     'TargetCategories','independent', ...
            % %                     'DataFormat','CB');
            %                 addpoints(lineLossValidation,iteration,double(gather(extractdata(lossValidation))))
            %
            %                 % Labeling F-score.
            %                   YPredValidation = extractdata(dlYPred_concat_test) == max(dlYPred_concat_test);
            %                   score = labelingFScore(YPredValidation,Y_Test);
            %                   addpoints(lineFScoreValidation,iteration,extractdata(score))
            % %                 YPredValidation = extractdata(dlYPredValidation) > labelThreshold;
            % %                 score = labelingFScore(YPredValidation,TValidation);
            % %                 addpoints(lineFScoreValidation,iteration,double(gather(score)))
            % %
            %                 drawnow
            %
            %                 % przywracamy dropout, zeby dzialal podczas funkcji
            %                 % forward
            %                 lgraph = layerGraph(dlnet1);
            %
            %                 larray = [ dropoutLayer(0.1, 'Name','dropout_1') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_1',larray);
            %                 larray = [ dropoutLayer(0.1, 'Name','dropout_2') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_2',larray);
            %                 larray = [ dropoutLayer(0.1, 'Name','dropout_3') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_3',larray);
            %                 larray = [ dropoutLayer(0.1, 'Name','dropout_4') ];
            %                 lgraph = replaceLayer(lgraph,'dropout_4',larray);
            %
            %                 dlnet1 = dlnetwork(lgraph);
            %                 %}
            %
            %             end
            %         end
            %     end
        end
        %}
        
        %
        % model = mnrfit(features,label,'model','hierarchical');
        %% Save model
        
        
    end
end

save_12_ECG_model(dlnet1,output_directory,classes);

end
function save_12_ECG_model(model,output_directory,classes)
% Save results.
tmp_file = 'finalized_model_entry_official_1.mat';
filename=fullfile(output_directory,tmp_file);
save(filename,'model','classes','-v7.3');


disp('Done.')
end


% find unique number of classes
function classes = get_classes(input_directory,files)

classes={};
num_files = length(files);
k=1;
for i = 1:num_files
    g = strrep(files{i},'.mat','.hea');
    input_file = fullfile(input_directory, g);
    fid=fopen(input_file);
    tline = fgetl(fid);
    tlines = cell(0,1);
    
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
        if startsWith(tline,'#Dx')
            tmp = strsplit(tline,': ');
            tmp_c = strsplit(tmp{2},',');
            for j=1:length(tmp_c)
                idx2 = find(strcmp(classes,tmp_c{j}));
                if isempty(idx2)
                    classes{k}=tmp_c{j};
                    k=k+1;
                end
            end
            break
        end
    end
    
    fclose(fid);
    
end
classes=sort(classes);
end

function [data,tlines] = load_challenge_data(filename)

% Opening header file
fid=fopen([filename '.hea']);

if (fid<=0)
    disp(['error in opening file ' filename]);
end

tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

f=load([filename '.mat']);

try
    data = f.val;
catch ex
    rethrow(ex);
end

end


function fs = get_fs(hea_data)

fid=hea_data;
tline = fgetl(fid);
tmp_str = strsplit(tline,' ');

fs = str2num(tmp_str{3});

end


function [recording_label,classes_label,single_recording_labels, fs]=get_true_labels(input_file,classes)

classes_label=classes;
single_recording_labels=zeros(1,length(classes));

fid=fopen(input_file);
tline = fgetl(fid);
tmp_str = strsplit(tline,' ');
recording_label = tmp_str{1};

fs = str2num(tmp_str{3});

tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
    if startsWith(tline,'#Dx')
        tmp = strsplit(tline,': ');
        tmp_c = strsplit(tmp{2},',');
        for j=1:length(tmp_c)
            %                 	        idx2 = find(strcmp(classes,str2double(tmp_c{j})));
            idx2 = find(classes==str2double(tmp_c{j}));
            
            single_recording_labels(idx2)=1;
        end
        break
    end
end
fclose(fid);

end

%% HELPER FUNCTION
function [gradients1, loss,dlYPred_concat] = modelGradients_demo(dlnet1,dlX1,Y)

%     dlYPred1 = forward(dlnet1,dlX1);
%     dlYPred2 = forward(dlnet2,dlX2);
%     dlX_concat=[dlYPred1;dlYPred2];
%     dlX_concat=reshape(dlX_concat,[1 40, 1, 128]);%the value 128 corresponds the mini batch size
dlX_concat=dlarray(single(dlX1),'SSCB');
dlY_concat=forward(dlnet1,dlX_concat);
%     dlYPred_concat = softmax(dlY_concat);
dlYPred_concat = sigmoid(dlY_concat);
loss = crossentropy(dlYPred_concat,Y,'TargetCategories','independent');
[gradients1] = dlgradient(loss,dlnet1.Learnables);


%     dlYPred1 = forward(dlnet1,dlX1);
% %     dlYPred2 = forward(dlnet2,dlX2);
%     dlX_concat=dlYPred1;
%     dlX_concat=reshape(dlX_concat,[1 size(Y,1), 1, 128]);%the value 128 corresponds the mini batch size
%     dlX_concat=dlarray(single(dlX_concat),'SSCB');
%     dlY_concat=dlYPred1;
%     dlYPred_concat = softmax(dlY_concat);
%     loss = crossentropy(dlYPred_concat,Y);
% %     loss2 = array2table(loss);
% %     loss2.Variables
%     [gradients1] = dlgradient(mean(loss),dlnet1.Learnables);
end

function score = labelingFScore(Y,T)

numObservations = size(T,2);

scores = (2 * sum(Y .* T)) ./ sum(Y + T);
score = sum(scores) / numObservations;

end


function dxmappingscored = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DXMAPPINGSCORED = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  DXMAPPINGSCORED = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  dxmappingscored = importfile("C:\PhysioNetChallenge2020_repo\evaluation-2020-master\dx_mapping_scored.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 27-Jun-2020 15:33:42

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Dx", "SNOMEDCTCode", "Abbreviation", "CPSC", "CPSCExtra", "StPetersburg", "PTB", "PTBXL", "Georgia", "Total", "Notes"];
opts.VariableTypes = ["string", "double", "string", "double", "double", "double", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Dx", "Abbreviation", "Notes"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Dx", "Abbreviation", "Notes"], "EmptyFieldRule", "auto");

% Import the data
dxmappingscored = readtable(filename, opts);

end
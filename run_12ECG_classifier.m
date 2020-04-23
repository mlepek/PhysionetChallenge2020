function [score, label] = run_12ECG_classifier(data,header_data,classes, model)

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    score = ones([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    features = get_12ECG_features(data,header_data);
    
    % default model
    %score = mnrval(model,features);		
    %[~,idx] = max (score);
    %label(idx)=1;
    
    % model BioS 23-04-2020, 57 manual features, shallow nnet
    
    % changing Inf and -Inf to NaN
    features( features == Inf ) = NaN;
    features( features == -Inf ) = NaN;
    
    % normalizing features
    features = (features - model.meanToNorm)./model.stdToNorm;
    
    % if there are NaNs in features vector then replace these NaNs by mean
    % values of the corresponding features
    if sum( isnan(features) ) > 0
        features( isnan(features) ) = model.mean_features( isnan(features) ); 
    end
    
    % prediction is performed here, features must be a column vector
    raw_predictions = model.net( features(:) );
    raw_predictions = raw_predictions';
    
    raw_predictions( raw_predictions < 0.0 ) = 0.0;
    
    % here, I create labels by taking the maximal value of prediction
    labels_predicted = zeros( 1, length(raw_predictions) );
    [~, max_ind] = max(raw_predictions);
    labels_predicted(max_ind) = 1;
    
    % here, I create labels by thresholding, so the illness is found if the
    % probability is over the threshold
    labels_predicted_by_threshold = ( raw_predictions > 0.3 );
    
    % here, I combine the two above label approaches (OR), so: always, we have at
    % least one illness detected, sometimes more than one.
    labels_predicted = labels_predicted | labels_predicted_by_threshold;
    
    % Finally, we put our results to the objects provided by organizers
    score = raw_predictions;
    label = labels_predicted;
    
end




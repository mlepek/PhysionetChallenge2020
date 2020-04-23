function model = load_12ECG_model()

        filename='finalized_model_BioS_23-04-2020.mat';
        A=load(filename);
        %model=A.model;
        
        % loading BioS model 18-04-2020
        
        model.net = A.net;
        model.meanToNorm = A.meanToNorm;
        model.stdToNorm = A.stdToNorm;
        model.mean_features = A.mean_features;

end



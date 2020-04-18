% This script runs the "driver" script in a handy way

%data_dir = 'D:\PnC2020\PhysioNetChallenge2020_Training_CPSC\Training_WFDB';
%output_dir = 'D:\PnC2020\PhysioNetChallenge2020_Training_CPSC\output';

data_dir = 'D:\PnC2020\PnC2020_TrainingData\Training_WFDB';
output_dir = 'D:\PnC2020\PnC2020_TrainingData\output_18-04-2020';

driver(data_dir, output_dir)


evaluate_12ECG_score(data_dir, output_dir, 'A_final_model_18-04-2020.csv');
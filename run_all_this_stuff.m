
% ten skrypt pozwoli wygodnie uruchamiac train model i driver

training_data = 'D:\PnC2020\Official_Data';
model_directory = 'D:\Repozytoria\PhysionetChallenge2020\output';
test_data = 'D:\PnC2020\Official_Data';
test_outputs = 'D:\PnC2020\Official_Outputs';

train_model(training_data, model_directory);

driver(model_directory, test_data, test_outputs);

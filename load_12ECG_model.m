function model = load_12ECG_model(model_directory)
	
	out_file='bios_model_entry_1.mat';
	filename=fullfile(model_directory,out_file);
        A=load(filename);
        model=A;
	
end



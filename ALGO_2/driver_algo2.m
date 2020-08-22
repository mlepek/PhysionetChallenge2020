function driver_algo2(input_directory)
	% Find files.
    input_files = {};
    for f = dir(input_directory)'
        if exist(fullfile(input_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'mat')
            input_files{end + 1} = f.name;
        end
    end

    % Iterate over files.
    num_files = length(input_files);
    for i = 1:num_files
        disp(['    ', num2str(i), '/', num2str(num_files), '...'])
        
        % Load data.
        file_tmp=strsplit(input_files{i},'.');
        tmp_input_file = fullfile(input_directory, file_tmp{1});
        [data,header_data] = load_challenge_data(tmp_input_file);
        
        features = algo2(data, header_data);
    end

    disp('Done.')
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


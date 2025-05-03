function extract_and_save_vars(sourceFolder, destinationFolder)
    % Ensure destination folder exists
    if ~exist(destinationFolder, 'dir')
        mkdir(destinationFolder);
    end

    % Get a list of all .mat files in the source folder
    matFiles = dir(fullfile(sourceFolder, '*.mat'));

    % Loop through each file
    for k = 1:length(matFiles)
        % Full path to the current file
        sourceFile = fullfile(sourceFolder, matFiles(k).name);

        % Load only the specified variables
        try
            data = load(sourceFile, 'arcsC3', 'nodesC2', 'Data','vessel_details', 'TaperVes','ploton','Scale','maxDaughters');
        catch ME
            warning('Could not load variables from %s: %s', matFiles(k).name, ME.message);
            continue;
        end

        % Create new file name with 'Network_' prefix
        newFileName = ['Network_' matFiles(k).name];
        destFile = fullfile(destinationFolder, newFileName);

        % Save the selected variables
        try
            save(destFile, '-struct', 'data');
        catch ME
            warning('Could not save variables to %s: %s', destFile, ME.message);
        end
    end

    fprintf('Processing complete. Processed %d files.\n', length(matFiles));
end

function addPath(folderPath)
    % Check if the folder exists
    if ~isfolder(folderPath)
        error('The specified folder does not exist.');
    end
    
    % Add the folder and its subfolders to the MATLAB path
    addpath(genpath(folderPath));
    
    disp(['The folder and all its subfolders have been added to the path: ', folderPath]);
end
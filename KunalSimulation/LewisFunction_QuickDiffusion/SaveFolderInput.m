%Function to receive input for the folder to save figures to during
%simulation, validate the input, and create folder if necessary

function save_D = SaveFolderInput(D)

    %Ask user to enter folder name in path D to save figures
    folder  = input('Enter name for folder to save figures:\n', 's');
    save_D = strcat( D, folder);
    
    %Check if folder exists. If it does, ask the user if they want to
    %overwrite files in it or enter a new name
    while isfolder( save_D)
        folder_response = input('Folder already exists. Type y to use this folder, or enter a new name:\n', 's');
        if folder_response ~= 'y'
            folder = folder_response;
            save_D = strcat( D, folder);
        else
            break
        end   
    end

    %Create folder if it doesn't exist
    if ~isfolder(save_D)
        mkdir( save_D)
    end
end

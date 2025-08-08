clc; close all; clear;
%%
name = 'QinTA';
mod = [ name '/Table'];
folder_path = ['./' mod '/'];

Tab = table();
% Get a list of all .mat files in the folder
file_list = dir(fullfile(folder_path, '*.mat'));
%%
 for i = 1:length(file_list)
    
    % Get the file name
    file_name = file_list(i).name;
    
    % Load the contents of the file
    load(fullfile(folder_path, file_name));
    
    Tab = [Tab; Data];
   

 end

Tab = rmmissing(Tab);

save(['Tablas/' name], 'Tab', '-v7.3');
writetable(Tab,['Tablas/' name '.csv']);
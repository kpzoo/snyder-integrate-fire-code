% Use sortSSAfols to sort all the data folders of SSA data
clearvars
clc
close all

% Home and main folders to consider
thisDir = cd;
% Main folders to sort into SSAfol names
cd('ssa data');
repSSA = dir('ssa*');
nIndep = length(repSSA); % no. repeats
cd(thisDir);

% SSA data in folders gamma[i], i refers to gamma value
SSAfol = {'gamma5', 'gamma10', 'gamma20', 'gamma30'};
lenSSAfols = length(SSAfol);

% Loop through and do the sorting
for i = 1:nIndep
    % Folders listed not in numerical order so use ssa_[num] structure
    nameID = repSSA(i).name;
    findID = regexp(nameID, '_');
    id1 = nameID(findID+1:end); % numerical part after '_'
    
    for j = 1:lenSSAfols
        % Identifiers within each SSAfol setting
        identif = ['IPP_' id1 '_' num2str(j) '*'];
        % Do the sorting
        SSApath = ['ssa data/' repSSA(i).name];
        sortSSAfols(SSAfol{j}, SSApath, identif);
    end
    % Progress
    disp(['Processed: ' repSSA(i).name]);
end
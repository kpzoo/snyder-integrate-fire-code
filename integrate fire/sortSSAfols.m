% Function that sorts the SSA data into SSAfols
function sortSSAfols(SSAfol, SSApath, identif)

% Home
thisDir = cd;

% Create a folder structure for the ssa data defined by SSAfol
cd(SSApath);
% Get relevant data in folder using identifier
datFile = dir(identif);
if exist(SSAfol, 'dir') ~= 7
    % Folder does not exist unless 7 returned, make folder
    mkdir(SSAfol);
end
% Move all desired files to SSAfol
for i = 1:length(datFile)
    movefile(datFile(i).name, SSAfol);
end
cd(thisDir);
% Obtain run segments and format QB data into files of fixed QB number
function nruns = formatQBdata(runsIn, runsdiv, folder)

% Assumptions and modifications
% - loadname and savename are subfolders of mainFol
% - assumes same tbio for all simulations
% - opens an extra file if 2 files cannot compose runsdiv
% - reprocesses same file if runsdiv < runs in each file

% Home directory
thisDir = cd;

% Obtain run segments and re-organise QB data
if runsIn < runsdiv
    error('Too few desired runs');
else
    % Get file names to divide
    cd(folder.QBundiv);
    QBfiles = dir('*.mat');
    flen = length(QBfiles);
    cd(thisDir);
    
    % Loop through files and re-organise data into segments of size runsdiv
    % by combining data in files and save in QBsave folder
    nsave = runsIn/runsdiv;
    if nsave ~= round(nsave)
        error('The ratio is runs/runsdiv is not integer');
    end

    % Variables to count loop iterations, QBs and files
    nruns = 0;
    iter = 0;
    file = 1;
    remRuns = 0;
    
    while nsave > 0
        % Load data and save in runsdiv components
        cd(folder.QBundiv);
        load(QBfiles(file).name);
        cd(thisDir);
        
        % Determine how many QBs remain in file
        Iset = QBout.Iset1(remRuns+1:size(QBout.Iset1, 1), 1:size(QBout.Iset1, 2));
        neff = QBout.neff - remRuns;
        % Load the next file if there are not enough QBs to satisfy a
        % division or if too many QBs in file then save another file
        if neff < runsdiv
            % Obtain remaining QBs needed to complete a file to be saved
            remRuns = runsdiv - neff;
            clear QBout
            cd(folder.QBundiv);
            % Update file for next run and account for if none available
            file = file + 1;
            if file > flen
                assignin('base', 'fileErr', file);
                error('Not enough files to meet required runs');
            end
            load(QBfiles(file).name);
            cd(thisDir);
            
            % Get extra QBs and check made required number else open a
            % further file
            nAvail = size(QBout.Iset1, 1);
            if remRuns < nAvail
                % Sufficient QBs in file to make remRuns
                Iset = [Iset; QBout.Iset1(1:remRuns, 1:size(QBout.Iset1, 2))];
            else
                % Not enough QBs so open a further file
                remRuns2 = remRuns - nAvail;
                disp(['Remaining runs 2 = ' num2str(remRuns2)]);
                Iset = [Iset; QBout.Iset1(1:nAvail, 1:size(QBout.Iset1, 2))];
                % Increment file index to account for etra file opening
                file = file + 1;
                if file > flen
                    assignin('base', 'fileErr', file);
                    error('Not enough files to meet required runs');
                end
                cd(folder.QBundiv);
                load(QBfiles(file).name);
                cd(thisDir);
                % Obtain Iset composed from 3 files and update remRuns 
                Iset = [Iset; QBout.Iset1(1:remRuns2, 1:size(QBout.Iset1, 2))];
                remRuns = remRuns2;
            end
            if size(Iset, 1) ~= runsdiv
                assignin('base', 'IsetErr', Iset);
                assignin('base', 'fileErr', file);
                error('Incorrect addition of extra QBs from subsequent file');
            end
        else
            % If extra runs are in the current file than required then
            % extract these and open same file again
            remRuns = neff - runsdiv;
            clear QBout
            cd(folder.QBundiv);
            if file > flen
                assignin('base', 'fileErr', file);
                error('Not enough files to meet required runs');
            end
            load(QBfiles(file).name);
            cd(thisDir);
            Iset = Iset(1:runsdiv, 1:size(QBout.Iset1, 2));
            remRuns = runsdiv;
        end
        
        % Save the reorganised data with runsdiv runs, assumes that all
        % simulations have same tbio 
        iter = iter + 1;
        cd(folder.QBload);
        tbio = QBout.tbio;
        name = ['QBdiv' num2str(iter)];
        save(name, 'Iset', 'tbio');
        cd(thisDir);
        
        % Account for save file and count the new QBs added
        nruns = nruns + size(Iset, 1);
        nsave = nsave - 1;
        clear QBout
    end    
    % Check that nruns meets desired runs (note runsIn is used vs runs as
    % this would be overwritten when QBdata is loaded)
    if nruns ~= runsIn
        assignin('base', 'nrunsErr', nruns);
        assignin('base', 'runsErr', runs);
        error('The expected number of QBs is not achieved');
    end
end
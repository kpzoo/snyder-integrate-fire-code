% Simple function to obtain time of code in hrs mins sec
function normaliseTime(runtime)

if runtime < 60
    disp(['Run time = ' num2str(runtime) ' secs']);
end
if runtime > 60 && runtime < 3600
    disp(['Run time = ' num2str(runtime/60) ' mins']);
end
if runtime > 3600
    disp(['Run time = ' num2str(runtime/3600) ' hrs']);
end
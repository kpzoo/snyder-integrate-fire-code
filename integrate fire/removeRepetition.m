% Algorithm that slightly translates coincident photons to ensure the
% estimated photon stream has no repeat times
function [TnewEst, normVal, newlen, oldlen] = removeRepetition(TphEst)

% Separate the unique elements from the repetitons
nph = length(TphEst);
[Tuni, iduni , ~] = unique(TphEst);
idrep = setdiff(1:nph, iduni);
idrep = idrep';
Trep = TphEst(idrep);

% Add the smallest possible translation to the non-unique elements
eps = min(TphEst)/10000000;
Trep = Trep + eps;
Tnew = [Tuni; Trep];
TnewEst = sort(Tnew);

% Obtain new unique length and norm and display results
oldlen = length(Tuni);
newlen = length(unique(TnewEst));
disp(['Unique no. photons times altered from ' num2str(oldlen) ' to ' num2str(newlen)]);
normVal = norm(TnewEst - TphEst);
disp(['The norm difference in photon times is ' num2str(normVal)]);
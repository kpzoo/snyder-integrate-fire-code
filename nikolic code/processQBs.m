% Function to process QB data and includes the experimental data for
% comparison - takes some code from gui_averageQBall 
function QBout = processQBs(t, Iset, runs, N, fcut, tstep, tMset, PLCmaxset, PLCtotset,...
    DAGmaxset, Ntrpmaxset, Camaxset, Cafreemaxset)

% Assumptions and modifications
% - filters QBs and compares them to experimental Hardie data
% - commented out many unneeded variables

% Initialise variables
Iset1 = zeros(1,N);  % filtered bumps 
Iset2 = zeros(1,N);  % filtered and shifted bumps 
Imaxset = 0; 
L = 0;
t12 = 0;
time1 = 0;
time2 = 0;
time3 = 0;
av_tM = 0;
av_peakPLC = 0;
av_PLCtot = 0;
av_peakDAG = 0;
av_peakNtrp = 0;
av_peakCa = 0;
av_peakCafree = 0;
w = 1; 
p = 1; 
s = 1; 
q = 1;
neff = runs;         % number of bumps "not failed" (Q.E. = neff/runs)

% Filter noisy QBs which are in Iset
kspan = 1 + uint16(1000/(3.1415*fcut*tstep));
If(1:kspan-1) = 0;
iss = 1:kspan-1;
ksfloat = single(kspan);
aa = 1-single(iss)./ksfloat;
cc = 2/ksfloat;

% Parameters obtained from comments in Drosophila_Phototransduction.m
Ibump = 2;
Ilat = 0.6;
tmin = 10;

% Main processing loop
disp('Started processing');

for m = 1:runs
    Inew = Iset(m,:);
    % If = Filter(Inew,t,fcut);
    for ns = kspan:N
        sumI = Inew(ns);
        for is = 1:kspan-1
            sumI = sumI + aa(is)*Inew(ns-is);
        end
        If(ns) = cc*sumI;
    end
    [Imax, k50] = max(-If);
    %Imaxall(m) = Imax;

    % Select the QBs only with Imax greater than some min value above the
    % noise, set by Ibump
    if Imax > Ibump
        lat = 0;
        t1 = 0;
        t2 = 0;
        for k=1:N
            if -If(k) > Ilat && lat == 0 && k > tmin/tstep
                lat =(k-1)*tstep;
            end
            if -If(k) >= Imax/2 && t1 == 0
                t1 = (k-1)*tstep;
            end
            if -If(k) <= Imax/2 && k>k50 && t2 == 0
                t2 = (k-1)*tstep;
            end
        end
        Iset1(w,:) = If;
        Imaxset(p) = Imax;
        L(s) = lat;
        time1(s) = t1-lat;
        tpeak = (k50-1)*tstep;
        time2(s) = tpeak - t1;
        time3(s) = t2-tpeak;
        t12(q) = (t1+t2)/2;
        av_tM = av_tM + tMset(w);
        av_peakPLC = av_peakPLC + PLCmaxset(w);
        av_PLCtot = av_PLCtot + PLCtotset(w);
        av_peakDAG = av_peakDAG + DAGmaxset(w);
        av_peakNtrp = av_peakNtrp + Ntrpmaxset(w);
        av_peakCa = av_peakCa + Camaxset(w);
        av_peakCafree = av_peakCafree + Cafreemaxset(w);
        w = w+1; 
        p = p+1; 
        s = s+1; 
        q = q+1;
    else
        neff = neff-1;
    end
    disp(['Processed run: ' num2str(m) ' of ' num2str(runs)]);
end

% Shifting and averaging for single QB across runs
T12 = 30;
Is = zeros(1,N);
summ1 = 0;
if neff > 0
    for m1 = 1:neff
        Is = zeros(1,N);
        I = Iset1(m1,:);        
        kin = round(t12(m1) - T12)/tstep;        
        if kin >= 0
            for k = 1:N-kin
                Is(k) = I(k+kin);
            end
        else
            kin = -kin;
            I(1:kin) = 0;
            for k = 1+kin:N
                Is(k) = I(k-kin);
            end
        end        
        Iset2(m1,:) = Is;        
        summ1=summ1+m1;        
    end    
    if summ1 > 1
        Q = mean(Iset2);
    else
        Q = Is;
    end    
else    
    Q = Is;
end
latencies = L;

% Experimental results for QBs from Hardie
bump = load('bump_data0.txt');
QQ = bump(:,1);
num2 = length(QQ); % 65;
TT = 0:1:(num2-1);
T = 0:tstep:(num2-1);
QB = spline(TT,QQ,T);
Qmax = max(-QB);
T1 = 0;
T2 = 0;
result = zeros(1,num2);

for k=1:length(T)
    if -QB(k) >= Qmax/2 && T1 == 0
        T1 = (k-1)*tstep;
    end
    if -QB(k) <= Qmax/2 && T1>0 && T2 == 0
        T2 = (k-1)*tstep;
    end
end

TQ12 = (T1+T2)/2;
Tend = (length(T)-1)*tstep;
Iav = Q;
t1 = 0;
t2 = 0;
[Imax,k50] = max(-Iav);

for k=1:length(t)
    if -Iav(k) >= Imax/2 && t1 == 0
        t1 = (k-1)*tstep;
    end
    if -Iav(k) <= Imax/2 && k>k50 && t2 == 0
        t2 = (k-1)*tstep;
    end
end
t12 = (t1+t2)/2;
istep = 1/tstep;
kin = round(t12 - TQ12)/tstep;
kend = 1 + Tend/tstep;

Is=zeros(1,kend);
if kin>=0
    for k=1:kend
        Is(k) = Iav(k+kin);
    end
else
    kin = -kin;
    Iav(1:kin) = 0;
    for k=1+kin:kend
        Is(k) = Iav(k-kin);
    end
end

for i2 = 1:num2
    ir = (i2-1)*istep + 1;
    result(i2) = Is(ir);  % reduce the number of points to exp
end

% Completion of processing and assigning of outputs
% av_Imax = mean(Imaxset);
% sd_Imax = sqrt(cov(Imaxset));
% av_L = mean(L);
% sd_L = sqrt(cov(L));
% av_t1 = mean(time1);
% sd_t1 = sqrt(cov(time1));
% av_t2 = mean(time2);
% sd_t2 = sqrt(cov(time2));
% av_t3 = mean(time3);
% sd_t3 = sqrt(cov(time3));
% av_tM = av_tM/runs;
% av_peakPLC = av_peakPLC/neff;
% av_PLCtot = av_PLCtot/neff;
% av_peakDAG = av_peakDAG/neff;
% av_peakNtrp = av_peakNtrp/neff;
% av_peakCa = av_peakCa/neff;
%%av_peakCafree = av_peakCafree/neff;

% Obtain peak current and latency to peak (which becomed a peak time when
% Tx2's first element is added via overlap function later
[Ipk, idpk] = min(Iset1,[], 2);
tpkLat = t(idpk);
tpkLat = tpkLat';
tLat = latencies;
disp('Processing complete');

% Assign output structures
QBout.Iset1 = Iset1;
% QBout.Iset2 = Iset2;
QBout.tbio = t;
QBout.tLat = tLat';
QBout.tpkLat = tpkLat;
QBout.Ipk = Ipk;
QBout.neff = neff;
QBout.runs = runs;

% Obtain some useful plots and assign some outputs to workspace
plotQBs(t, Iset1, T, TT, QQ, Is, result, latencies, neff);
% assignin('base', 'Iset1', Iset1);
% assignin('base', 'Iset2', Iset2);
% assignin('base', 'tLat', tLat);

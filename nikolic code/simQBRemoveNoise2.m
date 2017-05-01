% Main Nikolic model, no gui and the ability to simulate multiple photons
% Focus on setting noise flags to control cascade stochasticity
function QBout = simQBRemoveNoise2(runs, biosimname, flags, saveLoc)

% Assumptions and modifications
% - other parameters match default values of param_set from original code
% - removed comments related to activation variables
% - the Kcalx to Dpip2 params were commented in gui_single
% - there were commented section in gui_singleQB_multi

%% Set main simulation parameters

% Flags that control the noise settings
flagM = flags(1);
flagG = flags(2);
flagP = flags(3);
flagD = flags(4);
flagT = flags(5);
disp(['Flag settings: ' num2str(flags)]);

% Higher kMA gives longer tau_M lifetime, smaller Kcam makes it more sensitive
% to Ca influx - the aim is to have tau_M ~25ms, but latency time ~45ms

% Parameters of model (equivalent to contents of param_set) with
% descriprions and default values <------------ check with GUI for defaults

C = 2.4;            % Conversion factor in molecules/uM
ro = 0.25;          % Round onset is at 0.5+ro, so round(x-ro)=1 if x=0.9 (not used in stochastic mode)
% eff = 0.9;        % Collision efficiency

Kcalx = 0.2;     % Calcium concentration for 50% activity of pump in mM
vpkc = 0.06;      % Maximum rate of TRP phosphorylation by PKC (for Apkc = 1)
Kpkc1_50 = 100;     % DAG amount for 50% DAG-induced activity of PKC on TRP in molecules
Dpip2 = 6;     % Diffusion constant for PIP2 molecule in membrane in um^2/sec

bet = 1;            % Mean rate of channel closing (only for stochastic mode)
tDAGdelay = 12;     % tDAGdelay = tDAGdelay1
P1max = 1.0;        % P1max = P1max8;
Dgi = 1.2;          % Dgi = Dgi10;
Dga = 1.5;          % Dga = Dga11;

Mo = 1;             % Total number of isomerized M molecules
Gtot = 100;         % Total number of G protein molecules
Ptot = 100;         % Total number of PLC molecules
PIPo = 3e3;         % Total number of PIP2 molecules in dark
Atot = 70;          % Total number of Arrestin molecules
Ntot = 100;         % Total number of NINAC molecules
Ntrpo = 25;         % Total number of TRP channels
PKCtot = 100;       % Total number of PKC molecules
Kninac_max = 3;     % Maximum value for association constant NINAC-Arr2 in uM^(-1)
Kcam = 0.01;        % Concentration of Ca for 50% activity of calmodulin in mM
beta1 = 10;         % Exponential parameter for the reaction rate of NINAC + Arr -> NINAC*Arr

tGDP = 5;           % Time for GDP-GTP exchange in ms
t1 = 1;             % Delay time for Galpha-GTP protein activation in ms
tpi = 0.70;            % Time for reaction GPLC + PIP2-> DAG in dark in ms
Kpi_50 = 0.05;      % Conc of Ca for 50% activity of PLC in mM
beta2 = 3.5;         % Exponential parameter for rate of PLC + PIP2 -> PLC + DAG
tPdark = 100;       % Time constant for decay of GPLC in dark in ms
tDdark = 80;       % Time constant for decay of DAG in dark in ms
beta3 = 1.0;        % Exponential parameter for reaction rate of GPLC*+GAP->G+PLC
beta4 = 4.5;        % Exponential parameter for rate of DAG + DGK -> PA ... -> PIP2
K50_gap = 0.1;      % Concentration of calcium for 50% activity of GAP in mM
K50_dgk = 0.3;      % Concentration of calcium for 50% activity of DGK in mM
t2 = 20;            % Delay time for Ca-induced activation of GAP
t3 = 30;            % Delay time for Ca-induced activation of DGK
t4 = 30;            % Delay time for action of PKC on TRP
vph = 0.4;          % Rate of TRP relaxation in dark (dephosphorylation) in ms^(-1)
beta5 = 4;          % Exponential parameter for the reaction TRPp -> TRP

Dmetarh = 0.0;      % Diffusion const of active metarhodopsin um^2/s
alpha1 = 2e-3;
alpha2 = 2e-3;
alpha3 = 1.5e-3;
Dca = 220;          % Diffusion constant for Ca ion in microvillus in um^2/sec
Dmg = 200;          % Diffusion constant for Mg ion in microvillus in um^2/sec
Dna = 650;          % Diffusion constant for Na ion in microvillus in um^2/sec
Dk = 1000;          % Diffusion constant for K ion in microvillus in um^2/sec

Smv = 0.27;         % Surface of microvillus in um^2
Vmv = 4.24e-18;     % Volume of microvillus in L
Snk = 9.6e-4;       % Surface of microvillus neck in um^2
Lnk = 0.06;         % Length of microvillus neck in um

n = 4;              % Number of subunits in TRP channel (MWC model parameter)
Kcamtrp = 6;
Kpkc2_50 = 1.0;     % Calcium conc for 50% Ca-induced activity of PKC on TRP in mM

wca = 0.877;        % Fraction of total Permeability due to Ca
wmg = 0.101;        % Fraction of total Permeability due to Mg
wna = 0.011;        % Fraction of total Permeability due to Na
wk = 0.011;         % Fraction of total Permeability due to K

Vm = -0.070;        % Holding membrane potential in V
F = 96500;          % Faraday constant in C/mol

Ca_o = 1.5;         % Physiological Ca extracellular concentration in mM
Mg_o = 4;           % Physiological Mg extracellular concentration in mM
Na_o = 120;         % Physiological Na extracellular concentration in mM
K_o = 5.0;          % Physiological K extracellular concentration in mM

Ca_i = 1.6e-4;      % Physiological Ca intracellular concentration in mM
Mg_i = 3;           % Physiological Mg intracellular concentration in mM
Na_i = 8.0;         % Physiological Na intracellular concentration in mM
K_i = 140;          % Physiological K intracellular concentration in mM

% Set deterministic or stochastic mode
% flagM = 1;
% flagG = 1;
% flagP = 1;    % defaults for the case of all noise <---------------------
% flagD = 1;
% flagT = 1;

% Set simulation time parameters, cutoff and buffer <======================
time = 180;
tstep = 0.025; % usually 0.025
fcut = 100;
buffer = load('buffer50.txt');    % nonlinear buffer, calmodulin 0.5mM

% Parameters missed above but present in param_set
kMA = 0.005;
Yb_dark = 0.0000003;
Yb_max = 0.000007;
kninac_max = 3;
KR = 0.34;
KT = 0.0025;
Icalx_sat = 12.0;

% Obtain simulation time points (all in ms)
p = nextpow2(time/tstep);
Npoints = 2^p;
tup = (Npoints-1)*tstep;
t = 0:tstep:tup;
s = length(t);
itDAGdl = round(tDAGdelay/tstep);

% Variables to store molecular numbers and other simulation outputs with
% the "c" label indicating continous and "d" label indicating discrete

% M*
M = zeros(1,s);
M(1) = Mo;
MM = M;
Arr_free = zeros(1,s);
tau_M  = zeros(1,s);
Arr_free(1) = 1;
tau_M(1) = 1/(Arr_free(1)*kMA);

% G*
Gc1 = zeros(1,s);
Gc2 = zeros(1,s);
Gd = zeros(1,s);

% GPLC*
GPLCc1 = zeros(1,s);
GPLCc2 = zeros(1,s);
GPLCd = zeros(1,s);
PLCtot = 0;

% DAG
PIP2loss = zeros(1,s);
DAGc1 = zeros(1,s);
DAGc2 = zeros(1,s);
DAGd = zeros(1,s);
dDAGdt = zeros(1,s);

% TRP
CH = zeros(1,Ntrpo);    % array of channels (0 or 1) for stochastic mode
Ntrp = zeros(1,s);
Nact = zeros(1,s);
Nact(1) = Ntrpo;
Nactc1 = zeros(1,s);
Nactc2 = zeros(1,s);
%Yb = zeros(1,s);
Yb_old = Yb_dark;

% PLC, GAP, DGK, PKC
Aplc = zeros(1,s);
Agap = zeros(1,s);
Adgk = zeros(1,s);
Apkc = zeros(1,s);

% Ca2+, Mg2+, Na+, K+
Ca = zeros(1,s);
Ca(1) = Ca_i;
dCadt = zeros(1,s);
Cafree = zeros(1,s);
Cafree(1) = Ca_i;
dCafreedt = zeros(1,s);
Mg = zeros(1,s);
Mg(1) = Mg_i;
Na = zeros(1,s);
Na(1) = Na_i;
K = zeros(1,s);
K(1) = K_i;

% Itrp, Icalx, Iq
Itrp = zeros(1,s);
Icalx = zeros(1,s);
Ica = zeros(1,s);
Img = zeros(1,s);
Ina = zeros(1,s);
Ik = zeros(1,s);

% Dynamic threshold (= DAG amount needed for 1 open channel)
T = zeros(1,s);

% Further storage variables <-------------------------------------
Xg = zeros(1,s);
Xgap = zeros(1,s);
Xdgk = zeros(1,s);
Xpkc = zeros(1,s);

vcoll_PIP2 = zeros(1,s);
tau_PIP2 = zeros(1,s);
vPIP2 = zeros(1,s);

% Extra parameters not in Nikolic <========================================
multiP = 0;     % codes for multiple photons
cheat = 0;      % codes for a previous m file of params
% runs = 100;       % controls the number of repeat runs

% Just use a previous set of parameters
if cheat
    load('param_set.mat');
    for k = 1:77
        eval(param_set{k});
    end
end

% Alteration to allow for set of photons (this is in index form)
tphoton = [0 200];
photID = 1;
if max(tphoton) > tup - 100 && multiP
    error('Last photon too late for simulation');
end

% Declare variables to store data across runs
Iset = zeros(runs, s);
if multiP
    % In multiple photon case first isomerisation not necessarily at t = 0
    % unless specified in tphoton
    M(1) = 0;
    tMset = zeros(runs, length(tphoton));
else
    % First isomerisation at t = 0 as in Nikolic original code
    M(1) = Mo;
    tMset = zeros(runs, 1);
end
PLCmaxset = zeros(runs, 1);
PLCtotset = zeros(runs, 1);
DAGmaxset = zeros(runs, 1);
Ntrpmaxset = zeros(runs, 1);
Camaxset = zeros(runs, 1);
Cafreemaxset = zeros(runs, 1);

%% Run actual simulations of the transduction cascade dynamics

disp('Started simulation');
for m = 1:runs
    for i = 2:s

        % RHODOPSIN DYNAMICS
        Acam = Ca(i-1)/(Ca(i-1)+Kcam);
        Kninac = kninac_max*exp(-beta1*Acam);
        a = Kninac;
        b = 1+ Kninac*(Ntot-Atot)/C;
        c = - Atot/C;
        Afree = round(C*(-b+sqrt(b^2-4*a*c))/(2*a)-ro);
        Arr_free(i) = Afree;
        tau_M(i) = 1/(kMA*Afree);

        if flagM == 0
            dMdt = - kMA*Afree*MM(i-1);
            MM(i) = MM(i-1) + dMdt*tstep;
            M(i) = round(MM(i)-ro);
        else
            deltaM = - poisson(kMA*Afree*M(i-1)*tstep);
            M(i) = M(i-1) + deltaM;
        end

        % Multiple photon option: to simulate another photon at time t_photon
        if multiP && (photID <= length(tphoton))
            % Add extra M isomerisations and preset times tphoton
            if i == tphoton(photID)/tstep 
                M(i) = Mo;
                photID = photID + 1;
            end
        end

        % Reset test
        if M(i) < 1e-4
            M(i) = 0;
        end

        % CASCADE DYNAMICS
        vcoll_Gi = alpha1*((Dmetarh+Dgi)/Smv)*(Gtot - Gd(i-1) - GPLCd(i-1));        % rate of G-M* coll ms^(-1)
        vcoll_Ga = alpha2*(Dga/Smv)*(Ptot-GPLCd(i-1))/((1+sqrt(GPLCd(i-1)/pi))^2);  % rate of G*-PLC coll ms^(-1)
        vcoll_PI = alpha3*(Dpip2/Smv)*(PIPo - PIP2loss(i-1)); % -DAGd(i-1));        % rate of PIP2-GPLC* coll ms^(-1)
        vcoll_PIP2(i) = vcoll_PI;

        Aplc(i) = Cafree(i-1)/(Cafree(i-1)+Kpi_50);
        tPI = tpi*exp(beta3*Aplc(i)); %
        tau_PIP2(i) = tPI;

        vGi = vcoll_Gi/(1+vcoll_Gi*tGDP);   % rate of Galpha generation in ms^(-1)
        vGa = vcoll_Ga;                     % rate of GPLC generation in ms^(-1)
        vPI = vcoll_PI/(1+vcoll_PI*tPI);    % rate of DAG generation in ms^(-1)
        vPIP2(i) = vPI;

        Agap(i) = Agap(i-1) + Xgap(i-1)*tstep;
        Sgap = dCadt(i-1)*(K50_gap/(K50_gap+Ca(i-1))^2);
        Xgap(i) = Xgap(i-1) + tstep*(Sgap - Xgap(i-1))/t2;

        Adgk(i) = Adgk(i-1) + Xdgk(i-1)*tstep;
        Sdgk = dCadt(i-1)*(K50_dgk/(K50_dgk+Ca(i-1))^2);
        Xdgk(i) = Xdgk(i-1) + tstep*(Sdgk - Xdgk(i-1))/t3;

        % Consistency test
        if Agap(i) < 0
            Agap(i) = 0;
            %istim2 = i;
        end
        if Adgk(i) < 0
            Adgk(i) = 0;
            %istim3 = i;
        end

        tP = tPdark*exp(-beta2*Agap(i));
        tD = tDdark*exp(-beta4*Adgk(i));

        if flagG == 0
            Xg(i) = Xg(i-1) + tstep*(vGi*M(i-1)-Xg(i-1))/t1;
            dGdt1 = Xg(i);
            %   dGdt1 = vGi*c1(i-1);
            dGdt2 = vGa*Gd(i-1);
            Gc1(i) = Gc1(i-1) + dGdt1*tstep;
            Gc2(i) = Gc2(i-1) + dGdt2*tstep;
            % 2 quantised or continuous deterministic values
            Gd(i) = round(Gc1(i)-ro) - round(Gc2(i)-ro);
            %   Gd(i) = Gc1(i) - Gc2(i);
        else
            Xg(i) = Xg(i-1) + tstep*(vGi*M(i-1)-Xg(i-1))/t1;
            dGdt1 = Xg(i);
            %   dGdt1 = vGi*c1(i-1);
            dGdt2 = vGa*Gd(i-1);
            deltaG1 = poisson(dGdt1*tstep);
            deltaG2 = poisson(dGdt2*tstep);
            if deltaG2 > Gd(i-1)+deltaG1
                deltaG2 = Gd(i-1)+deltaG1;
            end
            Gd(i) = Gd(i-1) + deltaG1 - deltaG2;
        end

        % Consistency test
        if Gd(i) < 0
            Gd(i) = 0;
        end

        % Reset test
        if Gd(i)==0 && M(i)==0
            Gc1(i) = 0;
            Gc2(i) = 0;
        end

        if flagP == 0
            dPdt1 = dGdt2;
            dPdt2 = (1/tP)*GPLCd(i-1);
            GPLCc1(i) = GPLCc1(i-1) + dPdt1*tstep;
            GPLCc2(i) = GPLCc2(i-1) + dPdt2*tstep;
            % 3 quantised or continuous deterministic values (G-PLC)
            GPLCd(i) = round(GPLCc1(i)-ro) - round(GPLCc2(i)-ro);
            %   GPLCd(i) = GPLCc1(i) - GPLCc2(i);
            % new line for counting total active PLC
            PLCtot=PLCtot+(round(GPLCc1(i)-ro)-round(GPLCc1(i-1)-ro));
        else
            dPdt2 = (1/tP)*GPLCd(i-1);
            deltaP1 = deltaG2;
            deltaP2 = poisson(dPdt2*tstep);
            GPLCd(i) = GPLCd(i-1) + deltaP1 - deltaP2;
            % new line for counting total active PLC
            PLCtot=PLCtot+deltaP1;
        end

        % Consistency test
        if GPLCd(i) < 0
            GPLCd(i) = 0;
        end

        % Reset test
        if GPLCd(i) == 0 && M(i) == 0
            GPLCc1(i) = 0;
            GPLCc2(i) = 0;
        end

        if flagD == 0
            dDdt1 = vPI*GPLCd(i-1);
            dDdt2 = (1/tD)*DAGd(i-1);
            DAGc1(i) = DAGc1(i-1) + dDdt1*tstep;
            DAGc2(i) = DAGc2(i-1) + dDdt2*tstep;
            %old    DAGd(i) = round(DAGc1(i)-ro) - round(DAGc2(i-1)-ro);
            % 4 quantised or continuous deterministic values (DAG) ?????????????????
            DAGd(i) = round(DAGc1(i)-ro) - round(DAGc2(i)-ro);
            %    DAGd(i) = DAGc1(i) - DAGc2(i);
            dDAGdt(i) = (DAGd(i)-DAGd(i-1))/tstep;
        else
            dDdt1 = vPI*GPLCd(i-1);
            dDdt2 = (1/tD)*DAGd(i-1);
            DAGc1(i) = DAGc1(i-1) + poisson(dDdt1*tstep);
            DAGc2(i) = DAGc2(i-1) + poisson(dDdt2*tstep);
            %    DAGd(i) = round(DAGc1(i)-ro) - round(DAGc2(i-1)-ro);
            DAGd(i) = DAGc1(i) - DAGc2(i);
        end
        % new line for counting consumed PIP2
        %    PIP2loss(i) = PIP2loss(i-1) + DAGd(i);
        PIP2loss(i) = DAGc1(i);

        % Consistency test
        if DAGd(i) < 0
            DAGd(i) = 0;
        end

        % Reset test
        if DAGd(i) == 0 && M(i) == 0
            DAGc1(i) = 0;
            DAGc2(i) = 0;
        end

        % CHANNEL DYNAMICS
        L = 0;
        if i > itDAGdl
            L = DAGd(i-itDAGdl);
        end
        %L = DAGd(i-1);
        Yo = Yb_old;
        Atrp = ((1+KR*L).^n )./( (1+KR*L).^n + (1/Yo)*(1+KT*L).^n );
        alf = bet*Atrp/(1-Atrp);   % rate of channel opening

        if flagT == 0
            % 5 quantised or continuous deterministic values (Ntrp)
            Ntrp(i) = round(Nact(i-1)*Atrp-ro);
            %    Ntrp(i) = Nact(i-1)*Atrp;
        else
            alfatau = alf*tstep;
            betatau = bet*tstep;
            for j=1:Nact(i-1)
                if CH(j)==0
                    % CH(j) = poisson(alf*tstep);
                    % if CH(j)>0  % 1
                    %   CH(j) = 1;
                    % end
                    if alfatau > rand(1)  % 1
                        CH(j) = 1;
                    end
                else
                    % CH(j) = 1 - poisson(bet*tstep);
                    % if CH(j)<1  % 0
                    %   CH(j)=0;
                    % end
                    if betatau > rand(1)  % 1
                        CH(j) = 0;
                    end
                end
            end
            for j = Nact(i-1)+1:Ntrpo
                CH(j)=0;
            end
            Ntrp(i) = sum(CH);
        end

        Acam = Ca(i-1)/(Ca(i-1)+Kcam);
        % Positive feedback via TRP activation (new function)
        Acamtrp = Ca(i-1)/(Ca(i-1)+Kcamtrp);
        Yb_new = Yb_dark + (Yb_max-Yb_dark)*Acamtrp; % Acam
        vrelax = vph*exp(-beta5*Acam);

        if Apkc(i-1)>0
            dNdt1 = vpkc*PKCtot*Apkc(i-1)*Nact(i-1);  %  Ntrp(i-1)
        else
            dNdt1=0;
        end
        dNdt2 = vrelax*(Ntrpo - Nact(i-1));

        Nactc1(i) = Nactc1(i-1) + dNdt1*tstep;
        Nactc2(i) = Nactc2(i-1) + dNdt2*tstep;
        %Nactc1(i) = Nactc1(i-1) + poisson(dNdt1*tstep);
        %Nactc2(i) = Nactc2(i-1) + poisson(dNdt2*tstep);

        if flagT == 0
            Nact(i) = Ntrpo - round(Nactc1(i)-ro) + round(Nactc2(i)-ro);
            %    Nact(i) = Ntrpo - Nactc1(i) + Nactc2(i);
        else
            Nact(i) = Ntrpo - round(Nactc1(i)-ro) + round(Nactc2(i)-ro);
        end

        Apkc(i) = Apkc(i-1) + Xpkc(i-1)*tstep;
        A1 = (DAGd(i-1)/(Kpkc1_50 + DAGd(i-1)));
        A2 = (Cafree(i-1)/(Kpkc2_50 + Cafree(i-1)));
        Spkc = A1*dCafreedt(i-1)*Kpkc2_50/(Kpkc2_50+Cafree(i-1))^2 + A2*dDAGdt(i-1)*Kpkc1_50/(Kpkc1_50+DAGd(i-1))^2 ;
        Xpkc(i) = Xpkc(i-1) + tstep*(Spkc - Xpkc(i-1))/t4;

        %consistency test
        if Apkc(i) < 0
            Apkc(i) = 0;
            %istim4 = i;
        end

        % Threshold update
        alpha = (Nact(i)*Yb_new)^(1/4);
        if alpha==0
            T(i)=inf;
        else
            T(i) = (1-alpha)/(alpha*KR - KT);
        end

        %CURRENT DYNAMICS
        Icalx(i) = -Icalx_sat*Cafree(i-1)/(Cafree(i-1) + Kcalx);
        Ica(i) = ghk_current(Ntrp(i-1),Cafree(i-1),Ca_o,P1max,wca,2,Vm,Smv);
        Img(i) = ghk_current(Ntrp(i-1),Mg(i-1),Mg_o,P1max,wmg,2,Vm,Smv);
        Ina(i) = ghk_current(Ntrp(i-1),Na(i-1),Na_o,P1max,wna,1,Vm,Smv);
        Ik(i) = ghk_current(Ntrp(i-1),K(i-1),K_o,P1max,wk,1,Vm,Smv);

        Itrp(i) = Ica(i) + Img(i) + Ina(i) + Ik(i);
        Ica_leak = leak_current(2, Dca, Cafree(i-1), Ca_i, Lnk, Snk);
        Img_leak = leak_current(2, Dmg, Mg(i-1), Mg_i, Lnk, Snk);
        Ina_leak = leak_current(1, Dna, Na(i-1), Na_i, Lnk, Snk);
        Ik_leak = leak_current(1, Dk, K(i-1), K_i, Lnk, Snk);

        %ION DYNAMICS
        dCadt(i) = 1e-12*(-Ica(i-1)/2 + Icalx(i-1) - Ica_leak/2)/(F*Vmv);      % in mmol/(L*msec)
        Ca(i) = Ca(i-1) + dCadt(i)*tstep;                                      % in mM
        if Ca(i)<0
            Ca(i)=0;
        end
        i_buff = uint16(round(1e3*Ca(i)) + 1);
        B = 1 + buffer(i_buff); %/2;    % buffer power for 0.5 mM (0.25mM) Calmodulin)
        Cafree(i) = Ca(i)/B;          % in mM
        dCafreedt(i) = (Cafree(i)-Cafree(i-1))/tstep;
        dMgdt = 1e-12*(-Img(i-1) - Img_leak)/(2*F*Vmv);
        Mg(i) = Mg(i-1) + dMgdt*tstep;

        dNadt = 1e-12*(-Ina(i-1) - 3*Icalx(i-1) - Ina_leak)/(F*Vmv);
        Na(i) = Na(i-1) + dNadt*tstep;
        dKdt = 1e-12*(-Ik(i-1) - Ik_leak)/(F*Vmv);
        K(i) = K(i-1) + dKdt*tstep;

        % Consistency test
        if Ca(i) < 0
            Ca(i) = 0;
            Cafree(i) = 0;
        end
        Yb_old = Yb_new;

    end

    % Output values
    I = Itrp + Icalx;
    if ~multiP
        tM = sum(M)*tstep;
    else
        % Ensure tM gets all the deactivations of M
        idM = find(diff(M) == -1);
        tM = t(idM);
        if length(idM) ~= length(tphoton)
            error('Not all isomerisations recovered');
        end
    end
    PLCmax = max(GPLCd);
    Ntrpmax = max(Ntrp);
    Camax = max(Ca);
    Cafreemax = max(Cafree);
    DAGmax = max(DAGd);
    Imax = max(-I);

    % Output storage across runs
    Iset(m,:) = I;
    tMset(m,:) = tM;
    PLCmaxset(m) = PLCmax;
    PLCtotset(m) = PLCtot;
    DAGmaxset(m) = DAGmax;
    Ntrpmaxset(m) = Ntrpmax;
    Camaxset(m) = Camax;
    Cafreemaxset(m) = Cafreemax;
    disp(['Simulated run: ' num2str(m) ' of ' num2str(runs)]);
    
end
disp('All simulations complete');

%% Save and process data

% Post processing and plotting functions
N = s;
QBout = processQBs(t, Iset, runs, N, fcut, tstep, tMset, PLCmaxset, PLCtotset,...
    DAGmaxset, Ntrpmaxset, Camaxset, Cafreemaxset);
clearvars -except QBout biosimname totruns runs len biosimroot saveLoc flags

% Save processed data in correct sub-folder
thisDir = cd;
cd(saveLoc);
save([biosimname 'Pcs.mat'], 'biosimname', 'QBout', 'runs');
cd(thisDir);
% delete([biosimname 'Raw.mat']);
close all
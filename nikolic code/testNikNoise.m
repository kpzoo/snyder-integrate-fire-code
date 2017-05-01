% Examine the removal of different cascade noise components
clearvars
clc
close all

%% Setup code and run Nikolic model

% Basic test settings
runs = 1000;
biosimroot = 'QBnoise';

% Flags control noise [flagM flagG flagP flagD flagT]
flags = [0 0 0 1 1];
% Use this as identifier for saved file with no spaces
flagStr = num2str(flags);
flagStr = flagStr(~isspace(flagStr));

% Codes
% [0 0 0 0 0] -> shape and latency deterministic and singular

% Main Nikolic code
biosimname = [biosimroot flagStr];
QBout = simQBRemoveNoise(runs, biosimname, flags);

%% Extract and process outputs

% Get  basic outputs of currents, time and no. QBs
tbio = QBout.tbio;
Iset1 = QBout.Iset1;
neff = QBout.neff;

% Peak current and mean value
Ipk = QBout.Ipk;
meanPk = mean(Ipk);

% Current integral for QBs (across rows, hence 2)
QBarea = trapz(tbio, Iset1, 2);
meanArea = mean(QBarea);

% Latencies are shifted by 1ms as that's rhodopsin activation time 
tLat = QBout.tLat + 1;
meanLat = mean(tLat);


% Histogram of latencies
figure;
histogram(tLat, 20);
xlabel('latency (ms)');
ylabel('no. QBs');
box off
title(['Mean latency = ' num2str(meanLat) ' across ' num2str(neff) ' QBs']);

% Histogram of peak current
figure;
histogram(Ipk, 20);
xlabel('peak current (pA)');
ylabel('no. QBs');
box off
title(['Mean current = ' num2str(meanPk) ' across ' num2str(neff) ' QBs']);

% Histogram of area of QB
figure;
histogram(QBarea, 20);
xlabel('area (ms*pA)');
ylabel('no. QBs');
box off
title(['Mean area = ' num2str(meanArea) ' across ' num2str(neff) ' QBs']);


% Combine all plots
figure;
subplot(2, 2, 1);
histogram(tLat, 20);
xlabel('latency (ms)');
ylabel('no. QBs');
box off
title(['Mean latency = ' num2str(meanLat) ' across ' num2str(neff) ' QBs']);
subplot(2, 2, 2);
histogram(Ipk, 20);
xlabel('peak current (pA)');
ylabel('no. QBs');
box off
title(['Mean current = ' num2str(meanPk) ' across ' num2str(neff) ' QBs']);
subplot(2, 2, 3);
histogram(QBarea, 20);
xlabel('area (ms*pA)');
ylabel('no. QBs');
box off
title(['Mean area = ' num2str(meanArea) ' across ' num2str(neff) ' QBs']);
subplot(2, 2, 4);
plot(tbio, mean(Iset1), 'linewidth', 2);
box off
xlabel('time (ms)');
ylabel('current (pA)');
title(['Averaged raw QB (over ' num2str(neff) ' QBs)']);






% Basic average of QBs
figure;
plot(tbio, mean(Iset1), 'linewidth', 2);
box off
xlabel('time (ms)');
ylabel('current (pA)');
title(['Averaged raw QB (over ' num2str(neff) ' QBs)']);
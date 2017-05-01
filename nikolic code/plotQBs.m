% Function that plots QB data and includes the experimental data for
% comparison - takes some code from gui_averageQBall
function plotQBs(t, Iset1, T, TT, QQ, Is, result, latencies, neff)

% Plot shifted, averaged QB against experimental data
figure;
hold on
plot(TT,QQ,'or',T,Is,'k',TT,result,'x');
title(['Average QB (no.of bumps=',num2str(neff),') and experimental values (red circles)'])
xlabel('time [ms]')
ylabel('current [pA]')
hold off


% Latency distribution + 1ms (1ms estimated Rhodopsin activation time)
iendL = length(latencies);
xL1 = 0:10:100;
boxL1 = zeros(1,11);
r = 1;
for i=1:iendL
    if latencies(i)<5-r
        boxL1(1)=boxL1(1)+1;
        elseif latencies(i)<15-r
        boxL1(2)=boxL1(2)+1;
        elseif latencies(i)<25-r
        boxL1(3)=boxL1(3)+1;
        elseif latencies(i)<35-r
        boxL1(4)=boxL1(4)+1;
        elseif latencies(i)<45-r
        boxL1(5)=boxL1(5)+1;
        elseif latencies(i)<55-r
        boxL1(6)=boxL1(6)+1;
        elseif latencies(i)<65-r
        boxL1(7)=boxL1(7)+1;
        elseif latencies(i)<75-r
        boxL1(8)=boxL1(8)+1;
        elseif latencies(i)<85-r
        boxL1(9)=boxL1(9)+1;
        elseif latencies(i)<95-r
        boxL1(10)=boxL1(10)+1;
        elseif latencies(i)<105-r
        boxL1(11)=boxL1(11)+1;
    end
end
figure;
bar(xL1,boxL1,1,'g')
title('Distribution of latency times')
xlabel('time [ms]')
%ylabel('number')


% Basic figures for filtered QBs from each run
figure;
try
    plot(t', Iset1);
catch
    disp(size(t));
    disp(size(Iset1));
end
xlabel('t');
ylabel('filtered QBs');
title('Filtered QBs');

% Plot of sum of QBs (provided more than 1 effective run)
if neff > 1
    sumI = sum(Iset1);
    figure;
    plot(t, sumI);
    xlabel('t');
    ylabel('filtered QB sum');
    title('Sum of filtered QBs');
end

% Figure that includes the peak times
% Ipk = min(Iset1, [], 2);
% figure;
% plot(t, Iset1);
% hold on
% for i = 1:length(tMset)
%     plot([tMset(i) tMset(i)], [0 Ipk(i)]);
% end
% hold off
% xlabel('t');
% ylabel('filtered QBs');
% title('Filtered QBs with their maxima points');


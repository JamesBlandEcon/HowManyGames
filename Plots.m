%% Import the data
clear; clc; close all
rand('state',42)

set(0,'defaulttextInterpreter','latex') 
set(0,'defaultlegendInterpreter','latex') 
set(0,'defaultAxesTickLabelInterpreter','latex')

load MATLABdata
who

VList = who;

usecolor = 0; % turn to 1 to use colored bar charts

RPeriod = Period>0; % restrict analysis to these periods
%% Get id lists for groups, periods, and subjects

GroupList = unique(ugroupid);
SubjList  = unique(uid);
PeriodList = unique(Period);
TreatmentList = unique(Treatment);

%% Convert some veriables into some simpler names

A = 1-ActionAB;
B = ActionAB;
C = 1-ActionCD;
D = ActionCD;

ActComb(:,1) = A.*C; ActComb(:,2) = A.*D; ActComb(:,3) = B.*C; ActComb(:,4) = B.*D;

%% Choice probabilities by group

CPgroupB = zeros(numel(GroupList),4);
CPgroupN = zeros(numel(GroupList),2);

GroupTreatment = zeros(numel(GroupList),1);
for gg = 1:numel(GroupList)
    R = ugroupid==GroupList(gg) & RPeriod;
    CPgroupB(gg,:) = mean(ActComb(R,:));
    CPgroupN(gg,:) = mean([A(R) C(R)]) ;
    GroupTreatment(gg) = mean(Treatment(R));
end

CPBm = zeros(numel(TreatmentList),4);
CPNm = zeros(numel(TreatmentList),2);
CPBse = zeros(numel(TreatmentList),4);
CPNse = zeros(numel(TreatmentList),2);
for tt = 1:numel(TreatmentList)
    ntt = sum(GroupTreatment == TreatmentList(tt));
    CPBm(tt,:) = mean(CPgroupB(GroupTreatment == TreatmentList(tt),:));
    CPNm(tt,:) = mean(CPgroupN(GroupTreatment == TreatmentList(tt),:));
    CPBse(tt,:) = std(CPgroupB(GroupTreatment == TreatmentList(tt),:))./sqrt(ntt-1);
    CPNse(tt,:) = std(CPgroupN(GroupTreatment == TreatmentList(tt),:))./sqrt(ntt-1);
end


%% Compare actions: Broad analisis



h = figure;
hold all
if usecolor==0
    colormap(gray);
end
%errorbarbar(1:4,CPBm',CPBse');

hb = bar(1:4, CPBm');
pause(0.1);
for tt = 1:3
   xData = hb(tt).XData+hb(tt).XOffset;
   errorbar(xData,CPBm(tt,:),CPBse(tt,:),'k.')
end

%errorbarbar((1:4)',CPBm(1,:)',CPBse(1,:)')
legend('T1','T2','T3','Location','NorthWest')
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',['AC';'AD';'BC';'BD'])
xlabel('Action')
ylabel('Choice rate')
whitespace
hold off
saveas(h,'outputs/BroadActions.png')


%% Compare actions: Narrow analysis


h = figure;
hold all
if usecolor==0
    colormap(gray);
end
hb = bar(1:2, CPNm');
pause(0.1);
for tt = 1:3
   xData = hb(tt).XData+hb(tt).XOffset;
   errorbar(xData,CPNm(tt,:),CPNse(tt,:),'k.')
end
%errorbarbar(1:2,CPNm',CPNse');
legend('T1','T2','T3','Location','NorthWest')
set(gca,'XTick',[1:2])
set(gca,'XTickLabel',['A';'C'])
xlabel('Action')
ylabel('Choice rate')
whitespace
hold off
saveas(h,'outputs/NarrowActions.png')


%% By period

APeriod = zeros(numel(PeriodList),numel(TreatmentList));
CPeriod = zeros(numel(PeriodList),numel(TreatmentList));
for pp = 1:numel(PeriodList)
    for tt = 1:numel(TreatmentList)
        R = Period == PeriodList(pp) & Treatment == TreatmentList(tt); 
        APeriod(pp,tt) = mean(A(R));
        CPeriod(pp,tt) = mean(C(R));
    end
end

cmap = {'k','b','r'};

h = figure;
hold all
subplot(2,1,1)
hold all
for tt = 1:numel(TreatmentList)
    plot(PeriodList,APeriod(:,tt),strcat('-',cmap{tt}))
end
for tt = 1:numel(TreatmentList)
    plot([min(Period) max(Period)],mean(APeriod(:,tt))*[1;1],strcat('--',cmap{tt}))
end
hold off
ylim([0 1])
xlabel('Period')
ylabel('Pr[A]')
legend('T1','T2','T3','Location','best')
subplot(2,1,2)
hold all
for tt = 1:numel(TreatmentList)
    plot(PeriodList,CPeriod(:,tt),strcat('-',cmap{tt}))
end
for tt = 1:numel(TreatmentList)
    plot([min(Period) max(Period)],mean(CPeriod(:,tt))*[1;1],strcat('--',cmap{tt}))
end
hold off
ylim([0 1])
xlabel('Period')
ylabel('Pr[C]')
hold off

%% Compare predictions: Broad analisis



h = figure;
hold all
if usecolor==0
    colormap(gray);
end
%errorbarbar(1:4,CPBm',CPBse');

hb = bar(1:4, CPBm');
pause(0.1);
for tt = 1:3
   xData = hb(tt).XData+hb(tt).XOffset;
   errorbar(xData,CPBm(tt,:),CPBse(tt,:),'k.')
end

%errorbarbar((1:4)',CPBm(1,:)',CPBse(1,:)')
legend('T1','T2','T3','Location','NorthWest')
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',['AC';'AD';'BC';'BD'])
xlabel('Action')
ylabel('Choice rate')
whitespace
hold off
saveas(h,'outputs/BroadActions.png')



%% Compare predictions of learning model with 
clear; clc; clear global;

load C:\Users\jbland\Documents\Roommates\prPredictedLearning.mat

%%
set(0,'defaulttextInterpreter','latex') 
set(0,'defaultlegendInterpreter','latex') 
set(0,'defaultAxesTickLabelInterpreter','latex')

who

%%
ActionII = zeros([size(ACTCOMB) 4]);

for kk = 1:4
    ActionII(:,:,kk) = ACTCOMB==kk;
end

ActionII(:,:,5) = ActionII(:,:,1)+ActionII(:,:,2);
ActionII(:,:,6) = ActionII(:,:,1)+ActionII(:,:,3);
ActionII(:,:,7) = ActionII(:,:,2)+ActionII(:,:,3);

CPSIM(:,:,5,:) = CPSIM(:,:,1,:)+CPSIM(:,:,2,:);
CPSIM(:,:,6,:) = CPSIM(:,:,1,:)+CPSIM(:,:,3,:);
CPSIM(:,:,7,:) = CPSIM(:,:,2,:)+CPSIM(:,:,3,:);

CPSIM_B(:,:,5,:) = CPSIM_B(:,:,1,:)+CPSIM_B(:,:,2,:);
CPSIM_B(:,:,6,:) = CPSIM_B(:,:,1,:)+CPSIM_B(:,:,3,:);
CPSIM_B(:,:,7,:) = CPSIM_B(:,:,2,:)+CPSIM_B(:,:,3,:);
CPSIM_N(:,:,5,:) = CPSIM_N(:,:,1,:)+CPSIM_N(:,:,2,:);
CPSIM_N(:,:,6,:) = CPSIM_N(:,:,1,:)+CPSIM_N(:,:,3,:);
CPSIM_N(:,:,7,:) = CPSIM_N(:,:,2,:)+CPSIM_N(:,:,3,:);
%% pooled predictions

cspec = {'k','b','r','m','k','b'};

ActionLabels = {'$AC$','$AD$','$BC$','$BD$','$A$','$C$','$AD$ or $BC$'};

h = figure;
leg = {};
hold all
    for kk = 1:4
       plot(1:20,permute(mean(ActionII(:,:,kk),1),[3 2 1]),['--' cspec{kk}]);
         plt(kk) = plot(1:20,permute(mean(mean(CPSIM(:,:,kk,:),4),1),[3 2 1]),['-' cspec{kk}]);
         leg{numel(leg)+1} = ActionLabels{kk};
    end
    legend(plt,leg,'Location','best')
        xlabel('Period')
    ylabel('Choice rate')
    whitespace

hold off

h = figure;
leg = {};
hold all
    for kk = 5:6
       plot(1:20,permute(mean(ActionII(:,:,kk),1),[3 2 1]),['--' cspec{kk}]);
         plt(kk) = plot(1:20,permute(mean(mean(CPSIM(:,:,kk,:),4),1),[3 2 1]),['-' cspec{kk}]);
         leg{numel(leg)+1} = ActionLabels{kk};
    end
    legend(plt,leg,'Location','best')
        xlabel('Period')
    ylabel('Choice rate')
    whitespace

hold off

%%

UsePeriod = 1;

group = GROUP(:,1);
groupList = unique(group);

PREDICTION = zeros(numel(groupList),20-UsePeriod+1,numel(ActionLabels));
PREDICTION_B = zeros(numel(groupList),20-UsePeriod+1,numel(ActionLabels));
PREDICTION_N = zeros(numel(groupList),20-UsePeriod+1,numel(ActionLabels));
ACTUAL     = zeros(numel(groupList),20-UsePeriod+1,numel(ActionLabels));

for gg = 1:numel(groupList)
   PREDICTION(gg,:,:) = mean(mean(CPSIM(group==groupList(gg),UsePeriod:end,:,:),4),1);
   PREDICTION_B(gg,:,:) = mean(mean(CPSIM_B(group==groupList(gg),UsePeriod:end,:,:),4),1);
   PREDICTION_N(gg,:,:) = mean(mean(CPSIM_N(group==groupList(gg),UsePeriod:end,:,:),4),1);
   ACTUAL(gg,:,:)     = mean(ActionII(group==groupList(gg),UsePeriod:end,:),1);
end


% %%
% 
% for kk = 1:numel(ActionLabels)
% 
%     prediction = mean(PREDICTION(:,:,kk),2);
%     prediction_B = mean(PREDICTION_B(:,:,kk),2);
%     prediction_N = mean(PREDICTION_N(:,:,kk),2);
%     actual     = mean(ACTUAL(:,:,kk),2);
%     
%     b = regress(actual(:),prediction(:)); y = b*prediction(:);
%     b = regress(actual(:),prediction_B(:)); y_B = b*prediction_B(:);
%     b = regress(actual(:),prediction_N(:)); y_N = b*prediction_N(:);
% 
%     for bb = 1:numel(binsM)
%         meanChoice(bb) = mean(actual(prediction>=binsL(bb) & prediction<=binsR(bb)));
%     end
% 
%     h = figure;
%     m = max(max(prediction(:),actual(:)));
%     hold all
%         plot(prediction(:),actual(:),'.k')
%         plot([0 m],[0 m],'--k')
%         plot(binsM,meanChoice,'-xk')
%         plot(prediction(:),y,'-r')
%         xlabel(['Predicted choice frequency: ' ActionLabels{kk}])
%         ylabel(['Actual choice frequency:' ActionLabels{kk}])
%         whitespace
%     hold off
% end

%%
kmap = [1 2 4 5];


h = figure;
h.Position(3) = h.Position(3)*1.2
for kk = 1:4

    prediction = mean(PREDICTION(:,:,kk),2);
    prediction_B = mean(PREDICTION_B(:,:,kk),2);
    prediction_N = mean(PREDICTION_N(:,:,kk),2);
    actual     = mean(ACTUAL(:,:,kk),2);
    
    X = [ones(size(prediction(:))) log(prediction(:))-log(1-prediction(:))];
    [sX,II] = sort(X(:,2));
    [b,~,~,~,stats] = regress(log(actual(:))-log(1-actual(:)),X); yhat = X*b;
    y = 1./(1+exp(-yhat));
    r2 = stats(1);
        xylims = [min(min(prediction(:),actual(:))) max(max(prediction(:),actual(:)))];
    
    sp(kmap(kk)) =  subplot(2,3,kmap(kk))
    hold all
        plot(prediction(:),actual(:),'.k')
        %m = [min(prediction(:)) max(prediction(:))];
        plot(xylims,xylims,'--k')
        plot(prediction(II),y(II),'-k','LineWidth',2)
        %plot(prediction_B(:),y_B,'-r')
        %plot(prediction_N(:),y_N,'-b')
        if kmap(kk)==5
            xlabel(['Predicted choice rate (learning model)'])
        end
        %if kk == 1 | kk == 3
        %    ylabel(['Actual choice frequency'])
        %end
        title([ActionLabels{kk} ', $R^2=$' num2str(r2,'%1.2f')])
        xlim(xylims)
        ylim(xylims)
        whitespace
    hold off
end
subplot(2,3,3)
kk = 7
prediction = mean(PREDICTION(:,:,kk),2);
    prediction_B = mean(PREDICTION_B(:,:,kk),2);
    prediction_N = mean(PREDICTION_N(:,:,kk),2);
    actual     = mean(ACTUAL(:,:,kk),2);
    xylims = [min(min(prediction(:),actual(:))) max(max(prediction(:),actual(:)))];
     X = [ones(size(prediction(:))) log(prediction(:))-log(1-prediction(:))];
     [sX,II] = sort(X(:,2));
    [b,~,~,~,stats] = regress(log(actual(:))-log(1-actual(:)),X); yhat = X*b;
    y = 1./(1+exp(-yhat));
    r2 = stats(1);
    hold all
        plot(prediction(:),actual(:),'.k')
        plot(xylims,xylims,'--k')
        plot(prediction(II),y(II),'-k','LineWidth',2)
        %plot(prediction_B(:),y_B,'-r')
        %plot(prediction_N(:),y_N,'-b')
        %if kk == 1 | kk == 3
        %    ylabel(['Actual choice frequency'])
        %end
        title([ActionLabels{kk} ', $R^2=$' num2str(r2,'%1.2f')])
        xlim(xylims)
        ylim(xylims)
        whitespace
    hold off
    plt = subplot(2,3,6);
    plt.XTick = [];
    plt.YTick = [];
    plt.XColor = [1 1 1];
    plt.YColor = [1 1 1];
    hold all
    axis off
        p1 = plot(-1,-1,'.k')
        p2 = plot(-1,-1,'-k','LineWidth',2)
        p3 = plot(-1,-1,'--k')
        legend([p1 p2 p3],{'Predictions \& data','OLS fit','$45^\circ$ line'},'Box','off')
        ylim([0 1])
        xlim([0 1])
        
        
    hold off
    
    p1=get(sp(1),'position');
p2=get(sp(4),'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Actual choice rate (within-group mean)','visible','on');
saveas(h,'outputs/ForecastAccuracy.png')


%%
h = figure;
for kk = 5:6
    pp = kk-4;

    prediction = mean(PREDICTION(:,:,kk),2);
    prediction_B = mean(PREDICTION_B(:,:,kk),2);
    prediction_N = mean(PREDICTION_N(:,:,kk),2);
    actual     = mean(ACTUAL(:,:,kk),2);
    
    b = regress(actual(:),[ones(size(prediction(:))) prediction(:)]); y = [ones(size(prediction(:))) prediction(:)]*b;
    b = regress(actual(:),[ones(size(prediction(:))) prediction_B(:)]); y_B = [ones(size(prediction(:))) prediction_B(:)]*b;
    b = regress(actual(:),[ones(size(prediction(:))) prediction_N(:)]); y_N = [ones(size(prediction(:))) prediction_N(:)]*b;

    
     subplot(2,1,pp)
    hold all
        plot(prediction(:),actual(:),'.k')
        m = [min(prediction(:)) max(prediction(:))];
        plot(m,m,':k')
        plot(prediction(:),y,'-k')
        plot(prediction_B(:),y_B,'-r')
        plot(prediction_N(:),y_N,'-b')
        xlabel(['Predicted choice rate: ' ActionLabels{kk}])
        ylabel(['Actual choice rate:' ActionLabels{kk}])
        whitespace
    hold off
end
 
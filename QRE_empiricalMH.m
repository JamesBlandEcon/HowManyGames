%% Import the data
clear; clc; close all
clear global
rand('state',42)

load MATLABdata
who
ActList = {'AC','AD','BC','BD'};

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

%% Justification for ignoring first 10 periods

pAtt = zeros(numel(unique(Period)),3);
pCtt = zeros(numel(unique(Period)),3);

for pp =1:max(Period)
    for tt = 1:3
        pAtt(pp,tt) = mean(ActComb(Period==pp & Treatment ==tt,1),1)+mean(ActComb(Period==pp,2),1);
        pCtt(pp,tt) = mean(ActComb(Period==pp & Treatment ==tt,1),1)+mean(ActComb(Period==pp,3),1);
    end
end

cspec = {'k','r','b'};

%  h = figure;
%  hold all
%      for tt = 1:3
%          plot(1:max(Period),pAtt(:,tt),['-' cspec{tt}])
%          plot(1:max(Period),pCtt(:,tt),['--' cspec{tt}])
%      end
%      legend('p(A)','p(C)')
%  hold off

PeriodMin = 11;

for tt = 1:3
    ActionCount(tt,:) = sum(ActComb(Period>=PeriodMin & Treatment==tt,:),1);
end




%% set up some globals

global GB G1 G2 EmpPr NT d1 d2 ActComb uid ITREATMENT
drop = Period<11;
TrimList = {'ActComb','Treatment','a','c','x','y','d1','d2','uid','ugroupid'}
for vv = 1:numel(TrimList)
    eval([TrimList{vv} '(drop==1,:)=[];'])
end



EmpPr = zeros(size(ActComb));
%EmpPrGroup = ugroupid;
EmpPrGroup = Treatment;
PrGroupList = unique(EmpPrGroup);
for tt = 1:numel(PrGroupList)
    EmpPr(EmpPrGroup==PrGroupList(tt),:) = repmat(mean(ActComb(EmpPrGroup==PrGroupList(tt),:),1),[sum(EmpPrGroup==PrGroupList(tt)) 1]);
end

NT = size(EmpPr,1);




a = 10; c = 56; x = 100; y = 160;

GB = [a+c a+c a+c a+c;...
      a+y a   a+y a  ;...
      x+c x+c   c   c;...
      x+y x     y   0];
G1 = [a a a a;...
      a a a a  ;...
      x x 0 0;...
      x x 0 0];
G2 = [c c c c;...
      y 0 y 0;...
      c c c c;...
      y 0 y 0];

  % test
params = [0;0];
 llQREBroad(params)
 
 %% Simulation

 
 priorM = [0;-3;0]; priorV = eye(3);
 
 SimSize = 11000;
 BurnIn  = 1000;
 
 
 ModelList = {'Broad','Narrow','Mixture1','Mixture2'};
 nparams   = [2          2         3         3];
 
tt=1;
         ITREATMENT = Treatment==Treatment;
     
     for mm = 1:numel(ModelList)
          rand('state',42')
          randn('state',4242)
     eval(['llfun = @(x) llQRE' ModelList{mm} '(x)'])
     
     x = priorM;
     X = zeros(SimSize,numel(x));
     aLog = zeros(SimSize,1);
     likelihood = aLog;
     ax = llfun(x)+logmvnpdf(x',priorM',priorV);
     for ss = 1:SimSize

         
         prop = x+0.03*randn(size(x));
         aprop = llfun(prop)+logmvnpdf(prop',priorM',priorV);
         % note that since the proposal density is symmetric, we don't have to
         % worry about the other term
         log_a = aprop-ax; 
         accept = log_a>=log(rand(1,1));
         x = prop.*accept+x.*(1-accept);
         ax = aprop.*accept+ax.*(1-accept);
         X(ss,:) = x';
         aLog(ss) = accept;
         
         
         PriorDraw = mvnrnd(priorM,priorV)';
         
         likelihood(ss) = ax;

         disp([ModelList{mm} ' ' num2str(ss./SimSize)])

     end
     % evaluate marginal likelihood of model:
     % log(p(y|M)) = log(p(y|theta,M))+log(p(theta|M))-log(p(theta|y,M))
     %                  likelihood       prior for M      poseterior for M
     % since this does not change with theta (see Koop et al, Ex 16.7),
     % evaluate this at the posterior mode
     
     Xtrim = X((BurnIn+1):end,1:nparams(mm));
     postKdensity = mvksdensity(Xtrim,Xtrim,'Bandwidth',1.06.*std(Xtrim)./SimSize.^5);
     [F,II] = max(postKdensity);
     postModeX = Xtrim(II,:);
     lpyM(mm,tt) = llfun(postModeX')+logmvnpdf(postModeX,priorM(1:nparams(mm))',priorV(1:nparams(mm),1:nparams(mm)))-log(F);
     
     Sim(:,:,mm,tt)=X;
%      h = figure;
%      hold all
%         title([ModelList{mm} ' ' num2str(mean(aLog(BurnIn:end)))])
%         PlotThis = (X-repmat(mean(X),[SimSize 1]))./repmat(std(X),[SimSize 1]);
%         plot(PlotThis)
%         legend('r','lambda','mix')
%      hold off
     end

%%

save('C:\Users\jbland\Documents\Roommates\QRE.mat')
     
 
 %% Table
 clear 
 clear global
 load('C:\Users\jbland\Documents\Roommates\QRE.mat')
 
 Sim(1:BurnIn,:,:,:) = [];
 DATA = permute(Sim,[2 3 1 4]);

 
 PARAMS(1:2,:,:,:) =exp( DATA(1:2,:,:,:));
 PARAMS(3,:,:,:) = normcdf(DATA(3,:,:,:));
 fmt = '%0.3f';
 
 %%

 MNames = {'Broad','Narrow','Decisions','Subjects'};
  XLabels = {'$r$ - Risk preference','$\lambda$ - Choice precision','$\psi$ - Mixture'};
 fid = fopen(['outputs/QRETable.tex'],'w');

fprintf(fid,'\\begin{tabular}{l r@{.}l r@{.}l r@{.}l r@{.}l r@{.}l}\\\\ \\hline\\hline\n');
fprintf(fid,'&\\multicolumn{2}{c}{}&\\multicolumn{2}{c}{}&\\multicolumn{4}{c}{Mixture}\\\\ \n');
fprintf(fid,'Model');
for vv = 1:numel(ModelList)
    fprintf(fid,'&\\multicolumn{2}{c}{ %s }',MNames{vv});
end 

    fprintf(fid,' \\\\ \\hline \n');
for kk = 1:size(PARAMS,1)
    fprintf(fid,'%s',XLabels{kk});
    for vv = 1:numel(ModelList)
        if nparams(vv)>=kk
        v = mean(PARAMS(kk,vv,:,tt),3);
        integ=floor(abs(v));
        fract=abs(v)-integ;
        sign = '';
        if v<0
            sign = '-';
        end
        l = [sign int2str(integ)];
        %r = num2str(round(fract,decimals),decimals);
        r = num2str(fract,fmt);
        r = r(3:end);
        fprintf(fid,'& %s & %s',l,r);
        else
            fprintf(fid,'& \\multicolumn{2}{c}{}');
        end
        
    end
    fprintf(fid,' \\\\ \n');
    for vv = 1:numel(ModelList)
        if nparams(vv)>=kk
        v = std(permute(PARAMS(kk,vv,:,tt),[3 2 1]));
        integ=floor(abs(v));
        fract=abs(v)-integ;
        sign = '';
        if v<0
            sign = '-';
        end
        l = [sign int2str(integ)];
        %r = num2str(round(fract,decimals),decimals);
        r = num2str(fract,fmt);
        r = r(3:end);
        fprintf(fid,'& (%s & %s)',l,r);
        else
            fprintf(fid,'& \\multicolumn{2}{c}{}');
        end
    end
    fprintf(fid,' \\\\ \n');
end
fprintf(fid,' \\hline \n');
lpyR = lpyM(:,tt);
PostProb = round(exp(lpyR-max(lpyR))./sum(exp(lpyR-max(lpyR))),4);
fprintf(fid,'Posterior probability');
for mm = 1:numel(ModelList)
    v = PostProb(mm);
        integ=floor(abs(v));
        fract=abs(v)-integ;
        sign = '';
        if v<0
            sign = '-';
        end
        l = [sign int2str(integ)];
        %r = num2str(round(fract,decimals),decimals);
        r = num2str(fract,fmt);
        r = r(3:end);
        fprintf(fid,'& %s & %s',l,r);
end
fprintf(fid,' \\\\ ');

fprintf(fid,'\\hline\\hline\n');
fprintf(fid,' \\end{tabular}');
fclose(fid);


%%

rn = randn(10000,1);

h = figure;
hold all
    xlim([0.75 1.25])
    for mm = 1:numel(ModelList)
        [ff,xx] = ksdensity(exp(Sim(:,1,mm,1)));
        plt(mm) = plot(xx,ff);
    end
    [ff,xx] = ksdensity(exp(priorM(1)+rn.*sqrt(priorV(1,1))));
    plt(mm+1) = plot(xx,ff,'--k');
    MNames{mm+1} = 'prior';
    legend([plt],MNames)
    xlabel('r')
hold off

h = figure;
hold all
    xlim([0 0.1])
    for mm = 1:numel(ModelList)
        [ff,xx] = ksdensity(exp(Sim(:,2,mm,1)));
        plt(mm) = plot(xx,ff);
    end
    [ff,xx] = ksdensity(exp(priorM(2)+rn.*sqrt(priorV(2,2))));
    plt(mm+1) = plot(xx,ff,'--k');
    MNames{mm+1} = 'prior';
    legend([plt],MNames)
    xlabel('\lambda')
hold off

%% Predictions

mm = 2;

PrPredict = zeros(3,4,size(Sim,1));

for ss = 1:size(Sim,1)
    params = Sim(ss,:,mm);
    r = exp(params(1));
    l = exp(params(2));

    U = zeros(size(EmpPr));
    for aa = 1:4
        U(:,aa) = sum(EmpPr.*(repmat(G1(aa,:),[NT 1])+repmat(d1,[1 4])).^r,2)+sum(EmpPr.*(repmat(G2(aa,:),[NT 1])+repmat(d2,[1 4])).^r,2);
    end
    U = l.*(U-repmat(max(U,[],2),[1 4]));
    
    P = exp(U)./repmat(sum(exp(U),2),[1 4]);
    for tt = 1:3
        PrPredict(tt,:,ss) = mean(P(Treatment==tt,:));
    end
    disp(ss./size(Sim,1))
end 

%%

usecolor=0

load prPredictedLearning.mat

for tt = 1:3
    prEmpirical(tt,:) = mean(ActComb(Treatment==tt,:),1);
end

h = figure;
hold all
if usecolor==0
    colormap(gray);
end
%errorbarbar(1:4,CPBm',CPBse');


hb = bar(1:4, prEmpirical');
pause(0.1);
for tt = 1:3
   xData(tt,:) = hb(tt).XData+hb(tt).XOffset;
end
m = mean(PrPredict,3);

   pltQRE      = plot(xData(:),m(:),'sk','MarkerFaceColor',[1 1 1],'MarkerSize',9)
   pltLearning = plot(xData(:),prPredictedLearning(:),'dk','MarkerFaceColor',[1 1 1])


%errorbarbar((1:4)',CPBm(1,:)',CPBse(1,:)')
leg = legend([hb pltQRE pltLearning],{'T1 - empirical','T2 - empirical','T3 - empirical','QRE','Learning'},'Location','NorthWest','Box','off')
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',['AC';'AD';'BC';'BD'])
xlabel('Action')
ylabel('Choice rate (last 10 periods)')
hold off
saveas(h,'outputs/BroadActionsQRE.png')

%% Plot for learning model




%%

kill
%%
 h = figure;
hold all
if usecolor==0
    colormap(gray);
end
%errorbarbar(1:4,CPBm',CPBse');

prEmpiricalN = [prEmpirical(:,1)+prEmpirical(:,2) prEmpirical(:,1)+prEmpirical(:,3)];

hb = bar(1:2, prEmpiricalN');
pause(0.1);
xData = zeros(3,2)
for tt = 1:3
   xData(tt,:) = hb(tt).XData+hb(tt).XOffset;
end

m = mean(PrPredict,3);
mN = [m(:,1)+m(:,2) m(:,1)+m(:,3)];

prPredictedLearningN = [prPredictedLearning(:,1)+prPredictedLearning(:,2) prPredictedLearning(:,1)+prPredictedLearning(:,3)];
   pltQRE      = plot(xData(:),mN(:),'ob');
   pltLearning = plot(xData(:),prPredictedLearningN(:),'or');


%errorbarbar((1:4)',CPBm(1,:)',CPBse(1,:)')
leg = legend([hb pltQRE pltLearning],{'T1 - empirical','T2 - empirical','T3 - empirical','QRE','Learning'},'Location','best','Box','on')
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',['A';'C'])
xlabel('Action')
ylabel('Probability (last 10 periods)')
hold off
saveas(h,'outputs/NarrowActionsQRE.png')
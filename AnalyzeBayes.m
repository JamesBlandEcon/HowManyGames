%%
clear; clc;
modelname = 'ConstantOnly'
%modelname = 'Covariates'
filename = ['C:\Users\jbland\Documents\Roommates\MixModel' modelname '.mat']
eval(['load ' filename])
who *OUT

set(0,'defaulttextInterpreter','latex') 
set(0,'defaultlegendInterpreter','latex') 
set(0,'defaultAxesTickLabelInterpreter','latex')

%% some useful labels
%thetaLabels = {'r^{RD}','r^{LT}','\gamma','\kappa^{RD}','\kappa^{LT}'};
thetaLabels = {'r^{G}','r^{L}','\gamma','\kappa^{G}','\kappa^{L}'};
thetaTransform = {'log','log','normal cdf','probit','probit'};
taskLabels = {'G','L'};
if strcmp(modelname,'ConstantOnly')
    XLabels = {'Constant'};
end

%% remove burn in

BurnIn = 5000;

bOUT(:,:,1:BurnIn) = [];
lambdaOUT(1:BurnIn,:) = [];
wOUT(1:BurnIn) = [];
sigmaRRGOUT(:,:,1:BurnIn) = [];
thetaOUT(:,:,1:BurnIn) = [];
typePrOUT(:,:,1:BurnIn) = [];
typePriOUT(:,:,1:BurnIn) = [];

SimSize = size(lambdaOUT,1);

%% convert (transformed) theta to real theta
thetaOUT(:,1,:) = exp(thetaOUT(:,1,:));
thetaOUT(:,2,:) = exp(thetaOUT(:,2,:));
thetaOUT(:,3,:) = normcdf(thetaOUT(:,3,:));

%% marginal type probabilities

typePriM(:,:,1) = permute(sum(typePriOUT(:,3:4,:),2),[1 3 2]);
typePriM(:,:,2) = permute(sum(typePriOUT(:,[2 4],:),2),[1 3 2]);

typePrM(:,1) = mean(permute(sum(typePrOUT(:,3:4,:),2),[1 3 2]),1);
typePrM(:,2) = mean(permute(sum(typePrOUT(:,[2 4],:),2),[1 3 2]),1);

%%

a_sig = 0.05;

h = figure;
for tt = 1:3
    v = permute(thetaOUT(:,tt,:),[1 3 2]);
    mv = mean(v,2);
    qv = prctile(v,100.*[a_sig./2,1-a_sig./2],2);
    [temp,smv] = sort(mv);
    pv = abs(repmat(mv,[1 2])-qv);
    
    subplot(3,1,tt), errorbar((1:N)',mv(smv),pv(smv,1),pv(smv,2),'xk')
    hold on
    if tt<=2
        plot([0 N],[1 1],'--k')
        ylim([0 2])
    end
    xlim([0 N])
    ylabel(['$' thetaLabels{tt} '$'])
    whitespace
    hold off
end
saveas(h,['outputs/' modelname '_RRGDist.png'])

h = figure;
for tt = 1:2
    
    v = typePriM(:,:,tt);
    mv = mean(v,2);
    [temp,smv] = sort(mv);
    qv = prctile(v,100.*[a_sig./2,1-a_sig./2],2);
    [temp,smv] = sort(mv);
    pv = abs(repmat(mv,[1 2])-qv);
    
    PP{tt} = [mv pv];
    
    subplot(2,1,tt), errorbar((1:N)',mv(smv),pv(smv,1),pv(smv,2),'xk')
    xlim([0 N])
    ylabel(['$\Pr(\kappa^{' taskLabels{tt} '} = 1)$'])
end
saveas(h,['outputs/' modelname '_TypeDist.png'])

h = figure;
hold all

for ii = 1:N
    
    
    plot([PP{1}(ii,1),PP{1}(ii,1)],PP{2}(ii,1)+[-PP{2}(ii,2),PP{2}(ii,3)],'-r')
    plot(PP{1}(ii,1)+[-PP{1}(ii,2),PP{1}(ii,3)],[PP{2}(ii,1),PP{2}(ii,1)],'-r')
    
    
end
plot(PP{1}(:,1),PP{2}(:,1),'.k')
xlabel('$\Pr[\kappa^{RD}=1]$')
ylabel('$\Pr[\kappa^{LT}=1]$')
whitespace
hold off


%% construct mean and variances for transformed variables

mXB = zeros(1,5,SimSize);
sXB = zeros(1,3,SimSize);
for ss = 1:SimSize
    mXB(1,:,ss) = mean(X*bOUT(:,:,ss)',1);
    sXB(1,:,ss) = diag(sigmaRRGOUT(:,:,ss));
end

EXB = zeros(1,5,SimSize);
VXB = zeros(1,5,SimSize);



% first 2 are log normal
EXB(1,1:2,:) = exp(mXB(1,1:2,:)+sXB(1,1:2,:)./2);
VXB(1,1:2,:) = (exp(sXB(1,1:2,:))-1).*exp(2.*mXB(1,1:2,:)+sXB(1,1:2,:));

% gamma can be done with simulation
simG = 10000;
simGamma =normcdf(repmat(mXB(1,3,:),[simG 1 1])+repmat(sXB(1,3,:),[simG 1 1]).*repmat(randn(simG,1),[1 1 SimSize]));
EXB(1,3,:) = mean(simGamma,1);
VXB(1,3,:) = var(simGamma,[],1);

% last two are probits
EXB(1,4:5,:) = normcdf(mXB(1,4:5,:));
VXB(1,4:5,:) = normcdf(mXB(1,4:5,:)).*(1-normcdf(mXB(1,4:5,:)));


bOUT =permute(bOUT,[2 1 3]);
%% "Regression" output

reportstars = [0.05 0.01 0.001];

fmt = '%0.4f';

[M,SE,CI,STARS] = outputStats(bOUT,0.05,reportstars);

fid = fopen(['outputs/' modelname '_regOut.tex'],'w');

fprintf(fid,'\\begin{tabular}{l r@{.}l r@{.}l r@{.}l r@{.}l r@{.}l}\\\\ \\hline\\hline\n');
fprintf(fid,'Variable');
for vv = 1:numel(thetaLabels)
    fprintf(fid,'&\\multicolumn{2}{c}{$ %s $}',thetaLabels{vv});
end
fprintf(fid,' \\\\ \n');
fprintf(fid,'Transform');
for vv = 1:numel(thetaLabels)
    fprintf(fid,'&\\multicolumn{2}{c}{ %s }',thetaTransform{vv});
end
fprintf(fid,' \\\\ \\hline \n');
fprintf(fid,'\\multicolumn{11}{l}{{\\sc Mean coefficients $\\beta$}} \\\\  \n');
for kk = 1:size(bOUT,1)
    fprintf(fid,'%s',XLabels{kk});
    for vv = 1:numel(thetaLabels)
        v = M(kk,vv);
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
end
fprintf(fid,' \\\\ \n');
for kk = 1:size(bOUT,1)
    %fprintf(fid,'%s',XLabels{kk});
    for vv = 1:numel(thetaLabels)
        v = SE(kk,vv);
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
        fprintf(fid,'& (%s & %s)%s',l,r,STARS{kk,vv});
    end
end
fprintf(fid,' \\\\ \\hline \n');

fprintf(fid,'\\multicolumn{11}{l}{{\\sc Moments of transformed variables}');
if size(X,2)>1
    fprintf(fid,'(Evaluated at mean $X$)');
end
fprintf(fid,'} \\\\  \n');;
[M,SE,CI,STARS] = outputStats(EXB,0.05,reportstars);
    fprintf(fid,'Mean');
    for vv = 1:numel(thetaLabels)
        v = M(vv);
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
fprintf(fid,' \\\\ \n');
    %fprintf(fid,'%s',XLabels{kk});
    for vv = 1:numel(thetaLabels)
        v = SE(vv);
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
        fprintf(fid,'& (%s & %s)$^a$',l,r);
    end
fprintf(fid,' \\\\  \n');
[M,SE,CI,STARS] = outputStats(VXB,0.05,reportstars);
    fprintf(fid,'Variance');
    for vv = 1:numel(thetaLabels)
        v = M(vv);
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
fprintf(fid,' \\\\ \n');
    %fprintf(fid,'%s',XLabels{kk});
    for vv = 1:numel(thetaLabels)
        v = SE(vv);
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
        fprintf(fid,'& (%s & %s)$^a$',l,r);
    end
fprintf(fid,' \\\\ \\hline \n');

[M,SE,CI,STARS] = outputStats(sigmaRRGOUT,0.05,reportstars);

fprintf(fid,'\\multicolumn{11}{l}{{\\sc Covariance matrix} $\\Sigma$} \\\\  \n');
for rr = 1:3
    fprintf(fid,'$%s$',thetaLabels{rr});
    for cc = 1:3
        if rr>=cc
            v = M(rr,cc);
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
        elseif cc<3 
            fprintf(fid,'& \\multicolumn{2}{c}{ }');
        elseif rr<3
            fprintf(fid,'&\\multicolumn{2}{r}{$ %s $}',thetaLabels{rr+3});
            if rr ==1
                fprintf(fid,'&1&0000$^b$');
            elseif rr==2
                [wM,wSE,wCI,wSTARS] = outputStats(permute(wOUT,[3 2 1]),0.05,reportstars);
                v = wM;
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
                fprintf(fid,'&1&0000$^b$');
            end
        else
            fprintf(fid,'& \\multicolumn{2}{c}{ }');
        end
    end
    fprintf(fid,' \\\\  \n');
    for cc = 1:3
        if rr>=cc
            v = SE(rr,cc);
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
            star = STARS{rr,cc};
            if rr == cc
                star = '$^a$';
            end
            fprintf(fid,'& (%s & %s)%s',l,r,star);
        elseif cc<3 
            fprintf(fid,'& \\multicolumn{2}{c}{ }');
        elseif rr<3
            %fprintf(fid,'&\\multicolumn{2}{r}{$ %s $}',thetaLabels{rr+3});
            if rr ==1
                %
            elseif rr==2
                %[wM,wSE,wCI,wSTARS] = outputStats(permute(wOUT,[3 2 1]),0.05,reportstars);
                v = wSE;
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
            fprintf(fid,'&\\multicolumn{2}{c}{}& (%s & %s)%s',l,r,wSTARS{1});
           
                %fprintf(fid,'&1&0000$^b$');
            end
        else
            fprintf(fid,'& \\multicolumn{2}{c}{ }');
        end
    end
    fprintf(fid,' \\\\  \n');
end

fprintf(fid,' \\hline \n');
[M,SE,CI,STARS] = outputStats(permute(lambdaOUT,[2 3 1]),0.05,reportstars);
fprintf(fid,'\\multicolumn{11}{l}{{\\sc Logit choice precision}} \\\\  \n');
lLabels = {'$\\lambda^{G}$','$\\lambda^{L}$'};
for ll = 1:numel(lLabels)
    if ll==2
        fprintf(fid,'&\\multicolumn{2}{r}{');
    end
    fprintf(fid,lLabels{ll});
    if ll==2
        fprintf(fid,'}');
    end
    v = M(ll);
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
    %fprintf(fid,' & \n');
end
fprintf(fid,' \\\\ \n');
for ll = 1:numel(lLabels)
        if ll==2
            fprintf(fid,'&\\multicolumn{2}{r}{}');
        end
    v = SE(ll);
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
            fprintf(fid,'& (%s & %s)$^a$',l,r);
    
end
fprintf(fid,' \\\\ \n');

fprintf(fid,' \\hline\\hline \n');
% fprintf(fid,'\\multicolumn{11}{l}{\\footnotesize{');
% for ss = 1:numel(reportstars)
%     for tt = 1:ss
%     fprintf(fid,'*');
%     end
%     fprintf(fid,'$\\ p< %1.3f$,$\\quad$',reportstars(ss));
%     
% end
% fprintf(fid,'}}\\\\');
fprintf(fid,'\\multicolumn{11}{l}{\\footnotesize{$^*$ ($^{**}$, $^{***}$) A 95\\%% (99\\%%, 99.9\\%%) Bayesian credible region does not include zero}}\\\\');
fprintf(fid,'\\multicolumn{11}{l}{\\footnotesize{$^a$ Variable restricted to be positive. Stars supressed}}\\\\');
fprintf(fid,'\\multicolumn{11}{l}{\\footnotesize{$^b$ Variable restricted to one. Stars and standard deviation suppressed}}\\\\');
fprintf(fid,'\\end{tabular} \\\\ \n');
fclose(fid);

%% Correlation of r and testing p[broad|task]=p[broad]

corrR = reshape(sigmaRRGOUT(2,1,:)./sqrt(sigmaRRGOUT(1,1,:).*sigmaRRGOUT(2,2,:)),[],1,1);
disp('mean, std, p>0 of risk aversion correlation:')
disp([mean(corrR) std(corrR) mean(corrR>0)])

disp('posterior prob that broad bracketing more likely in RD:')
disp(mean(bOUT(:,4,:)>bOUT(:,5,:),3))
disp('95% CI of probability difference:')
disp(quantile(permute(normcdf(bOUT(:,4,:))-normcdf(bOUT(:,5,:)),[3,1,2]),[0.025,0.975],1))



%% Forgone expected utility

ExRate = 0.02;
FCE = zeros(N,20,SimSize);
FCEabs = zeros(N,20,SimSize);
safeCE = zeros(N,20,SimSize);
for ss = 1:SimSize
    [FEU,FCE(:,:,ss),FCEabs(:,:,ss),safeCE(:,:,ss)] = FEU_RD(thetaOUT(:,1,ss),thetaOUT(:,4,ss)>=0,thetaOUT(:,3,ss),lambdaOUT(ss,1),ACTCOMB,ACTCOMBPm,D1,D2);
    disp(ss)
end

% drop first 10 periods of calculations

dropPeriods = 10;
FCE(:,1:dropPeriods,:) = [];
FCEabs(:,1:dropPeriods,:) = [];
safeCE(:,1:dropPeriods,:) = [];
%%

mFCE = ExRate.*permute(mean(FCEabs,2),[1 3 2]);
PlotFCE = [mean(mFCE,2) prctile(mFCE,100*[a_sig./2 1-a_sig./2],2)];
[temp,sPFCE] = sort(PlotFCE(:,1));
h = figure;
hold all
    errorbar(1:N,PlotFCE(sPFCE,1),abs(PlotFCE(sPFCE,2)-PlotFCE(sPFCE,1)),abs(PlotFCE(sPFCE,3)-PlotFCE(sPFCE,1)),'xk')
    xlim([0 N])
    ylabel('Forgone  certainty equivalent (\$)')
    title(['mean value: $' num2str(mean(mFCE(:)),'%1.2f') '\ (' num2str(std(mFCE(:)),'%1.2f') ')$'])
    %ylim([0 1])
    whitespace
hold off
saveas(h,['outputs/' modelname 'CEpotential.png'])

mFCE = ExRate.*permute(mean(FCEabs.*repmat(thetaOUT(:,4,:)<0,[1 (20-dropPeriods) 1]),2),[1 3 2]);
PlotFCE = [mean(mFCE,2) prctile(mFCE,100.*[a_sig./2 1-a_sig./2],2)];
[temp,sPFCE] = sort(PlotFCE(:,1));
h = figure;
hold all
    errorbar(1:N,PlotFCE(sPFCE,1),abs(PlotFCE(sPFCE,2)-PlotFCE(sPFCE,1)),abs(PlotFCE(sPFCE,3)-PlotFCE(sPFCE,1)),'xk')
    xlim([0 N])
    title(['mean value: $' num2str(mean(mFCE(:)),'%1.2f') '\ (' num2str(std(mFCE(:)),'%1.2f') ')$'])
    %ylim([0 1])
    ylabel('Forgone  certainty equivalent (\$)')
    whitespace
hold off
saveas(h,['outputs/' modelname 'CEactual.png'])


%% Table of individual parameters

PrBG = sum(typePriOUT(:,3:4,:),2);
PrBL = sum(typePriOUT(:,[2 4],:),2);

FCEpotential = ExRate.*mean(FCEabs,2);
FCEactual = ExRate.*mean(FCEabs.*repmat(thetaOUT(:,4,:)<0,[1 (20-dropPeriods) 1]),2);


iParams = [thetaOUT(:,1:3,:) PrBG PrBL FCEpotential FCEactual];

iParamsMean = mean(iParams,3);
iParamsStd  = sqrt(var(iParams,[],3));

TreatmentN = Treatment(Period==20);

fid = fopen('outputs/iParamsTable.tex','w');

fprintf(fid,'\\begin{longtable}{rrrrrrrrr} \\hline\\hline \n');
fprintf(fid,' & &&  & & & \\multicolumn{2}{c}{Foregone CE (\\$)} \\\\ \n');
fprintf(fid,'Subject & $r^G$ & $r^L$ & $\\gamma$ & $\\rho^G$ & $\\rho^L$& Potential & Actual \\\\ \\hline');

for ii = 1:size(iParamsMean,1)
    fprintf(fid,int2str(ii));
    fprintf(fid,'&%1.3f',iParamsMean(ii,1:5));
    fprintf(fid,'&%1.2f',iParamsMean(ii,6:7));
    fprintf(fid,'\\\\');
    
    %fprintf(fid,int2str(ii));
    fprintf(fid,'&(%1.3f)',iParamsStd(ii,1:5));
    fprintf(fid,'&(%1.2f)',iParamsStd(ii,6:7));
    fprintf(fid,' \\\\ \\addlinespace[0.05cm] ');
    
end
fprintf(fid,' \\hline\\hline')

fprintf(fid,'\\end{longtable}');

fclose(fid);




%% Actual CE foregone and risk aversion

Rplot = [mean(thetaOUT(:,1,:),3)  permute(prctile(thetaOUT(:,1,:),100.*[a_sig./2 1-a_sig./2],3),[1 3 2])];

h = figure;
hold all
    plot(Rplot(:,1),PlotFCE(:,1),'.k')
    xlim([0 2])
    xlabel('r^G')
    ylabel('Certainty equivalent (fraction of optimal)')
    whitespace
hold off
saveas(h,['outputs/' modelname 'RvsCEactual.png'])

%% Predicted values of bar plots

CP = zeros(N,4,SimSize);

prBRD = typePriOUT(:,3,:)+ typePriOUT(:,4,:);

T(:,1) = D1(:,1)==0 & D2(:,1)==0;
T(:,2) = D1(:,1)==50 & D2(:,1)==0;
T(:,3) = D1(:,1)==0 & D2(:,1)==50;

CPSIM = zeros(N,20,4,SimSize);
CPSIM_B = zeros(N,20,4,SimSize);
CPSIM_N = zeros(N,20,4,SimSize);
for ss = 1:SimSize
    [cp,~,cpB,cpN] = RDChoiceProbs(thetaOUT(:,1,ss),prBRD(:,:,ss),thetaOUT(:,3,ss),lambdaOUT(ss,1),ACTCOMB,ACTCOMBPm,D1,D2);
    %for tt = 1:3     
    %    CP(tt,:,ss) = reshape(mean(mean(cp(T(:,tt),:,:),1),2),1,[],1);
    %end
    CP(:,:,ss) = permute(mean(cp(:,11:20,:),2),[1 3 2]);
    CPSIM(:,:,:,ss) = cp;
    CPSIM_B(:,:,:,ss) = cpB;
    CPSIM_N(:,:,:,ss) = cpN;
    disp(ss./SimSize)
end

%% Compare by subject

cspec = {'k','r','b'};
Predicted = zeros(3,20);
for aa = 1:4
    h = figure;
    hold all
    for tt = 1:3
        Predicted = mean(mean(CPSIM(T(:,tt),:,aa,:),4));
        Empirical = mean(ACTCOMB(T(:,tt),:)==aa,1);
        plot(1:20,Empirical,['--' cspec{tt}])
        plot(1:20,Predicted,['-' cspec{tt}])
    end
    
    hold off
end



%% compare individual actions
prEmpirical = zeros(3,4);
prPredicted = zeros(3,4);


for tt = 1:3
 for aa = 1:4
    prEmpirical(tt,aa) = mean(ActComb(Treatment==tt & Period>=11,1)==aa,1);
    prPredicted(tt,:) = mean(mean(CP(T(:,tt),:,:),3));
 end
 prPredicted(tt,:) = mean(mean(CP(T(:,tt),:,:),3));
end

CPLearning = CP;
TREATMENT = zeros(N,20);
GROUP = zeros(N,20);
for rr = 1:numel(uidList)
    for tt = 1:20
        TREATMENT(rr,tt) = Treatment(uid==uidList(rr) & Period==tt);
        GROUP(rr,tt) = ugroupid(uid==uidList(rr) & Period==tt);
    end
end


prPredictedLearning = prPredicted;
save('C:\Users\jbland\Documents\Roommates\prPredictedLearning.mat','prPredictedLearning','CPLearning','ACTCOMB','TREATMENT','GROUP','CPSIM','CPSIM_B','CPSIM_N')

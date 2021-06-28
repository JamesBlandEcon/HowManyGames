%% Mixture model 
clear; clc;

rand('state',42)
randn('state',4242)

load MATLABdata.mat

c1 = chi1(1,:);
c2 = chi2(1,:);


%% get LT data into right format

invest2 = (invest2_1+invest2_2)./2;

uidLT = uid(Period==20);
i1 = invest1(Period==20,:);
i2 = invest2(Period==20,:);

N = size(i1,1);

%% get RD data into the right format

ActComb = 1.*(ActionAB==0 & ActionCD==0)...
         +2.*(ActionAB==0 & ActionCD==1)...
         +3.*(ActionAB==1 & ActionCD==0)...
         +4.*(ActionAB==1 & ActionCD==1);


ACTCOMB = zeros(N,20); % identifier for the combination of actions taken by subject
D1 = zeros(N,20);
D2 = zeros(N,20);
ACTCOMBP = zeros(N,20,3); % Identifier for part

uidList = unique(uid);
ugroupidList = unique(ugroupid);

for ii = 1:N
    g = mean(ugroupid(uid==uidList(ii)));
    for pp = 1:20
        ACTCOMB(ii,pp)  = ActComb(uid == uidList(ii) & Period==pp);
        D1(ii,pp) = d1(uid == uidList(ii) & Period==pp);
        D2(ii,pp) = d2(uid == uidList(ii) & Period==pp);
        
        
        ACTCOMBP(ii,pp,:) = ActComb(ugroupid==g & uid~=uidList(ii) & Period==pp);
        
    end
end


ACTCOMBPm = zeros(N,20,4);
for aa = 1:4;
    ACTCOMBPm(:,:,aa) = mean(ACTCOMBP==aa,3);
end


%% generate interpolation grid for LT1 and LT2B

fR1 = @(x,c) log(c-1)./(log(c.*x+100-x)-log(100-x));

xgrid = linspace(1,99,99)';
rgrid1  = fR1(repmat(xgrid,[1 numel(c1)]),repmat(c1,[numel(xgrid),1]));

rgrid2N = fR1(repmat(xgrid,[1 numel(c2)]),repmat(c2,[numel(xgrid),1]));


%%

syms x r c

obj2 = @(x,r,c) (2.*c.*x + 200-2*x).^(1-r)./(1-r)...
            +2.*(1.*c.*x + 200-2*x).^(1-r)./(1-r)...
               +(200-2*x).^(1-r)./(1-r);

obj2sym = obj2(x,r,c);
obj2symFOC = diff(obj2sym,x);
obj2FOC = matlabFunction(obj2symFOC);

rgrid2B = zeros(size(rgrid2N));
for xx = 1:numel(xgrid)
    for cc = 1:numel(c2)
        obj = @(r) obj2FOC(c2(cc),r,xgrid(xx));
        rgrid2B(xx,cc) = fzero(obj,rgrid1(xx,cc));
    end
end

%% min distance estimates from LT (starting values)


rLT = zeros(N,1);
kLT = zeros(N,1);

for ii = 1:size(i1,1)
    objN = @(r) sum([i1(ii,:)-LTinterp(r,xgrid,rgrid1),i2(ii,:)-LTinterp(r,xgrid,rgrid2N),2].^2);
    [r1,o1] = fminsearch(objN,0.5);
    
    objB = @(r) sum([i1(ii,:)-LTinterp(r,xgrid,rgrid1),i2(ii,:)-LTinterp(r,xgrid,rgrid2B),2].^2);
    [r2,o2] = fminsearch(objB,0.5);
    
    if o2>o1
        rLT(ii) = r2;
        kLT(ii) = 1;
    else
        rLT(ii) = r1;
        kLT(ii) = 0;
    end
end

rLTminD =rLT;
kLTminD =kLT;




%% test LT predict function

LTpredict = LTinterpBoth(rLT,kLT,xgrid,rgrid1,rgrid2N,rgrid2B);

%% add a little bit to D1 and D2 so that earnings are always positive in the RD

D1 = D1; D2 = D2;

%% Min Distance estimator for kRD,rRD, gt
rRD = rLT;
kRD = kLT;

%for ii = 1:N
%    obj0 = @(r) -RDlike(r,0,0.8,0.05,ACTCOMB(ii,:),ACTCOMBPm(ii,:,:),D1(ii,:),D2(ii,:));
%    [r0,y0] = fminsearch(obj0,rRD(ii));
%    obj1 = @(r) -RDlike(r,1,0.8,0.05,ACTCOMB(ii,:),ACTCOMBPm(ii,:,:),D1(ii,:),D2(ii,:));
%    [r1,y1] = fminsearch(obj1,rRD(ii));
%    if y0>y1
%        rRD(ii) = r0;
%        kRD(ii) = 0;
%    else
%        rRD(ii) = r1;
%        kRD(ii) = 1;
%    end
%end

gt = zeros(N,1);
for ii = 1:N
    obj0 = @(x) -RDlike(x(1),0,x(2),0.05,ACTCOMB(ii,:),ACTCOMBPm(ii,:,:),D1(ii,:),D2(ii,:));
    [x0,y0] = fminsearch(obj0,[rRD(ii);log(0.2)]);
    obj1 = @(x) -RDlike(x(1),1,x(2),0.05,ACTCOMB(ii,:),ACTCOMBPm(ii,:,:),D1(ii,:),D2(ii,:));
    [x1,y1] = fminsearch(obj1,[rRD(ii);log(0.2)]);
    if y0>y1
        rRD(ii) = x0(1);
        gt(ii) = x0(2);
        kRD(ii) = 0;
    else
        rRD(ii) = x1(1);
        kRD(ii) = 1;
        gt(ii) = x1(2);
    end
    disp(ii)
end


RDmaxlike = [exp(rRD) normcdf(gt)];



%%
%typeFun = @(kRD,kLT) 0.*(kRD == 0 & kLT ==0) + 1.*(kRD == 0 & kLT ==1) + 2.*(kRD == 1 & kLT ==0) + 3.*(kRD == 1 & kLT ==1);
%type = typeFun(kRD,kLT);

%% Define characteristics vector X




X = ones(N,1);

%modelname = 'Covariates';
%modelname = 'CovariatesSmall'
%modelname = 'MFAC'
%modelname = 'Sex'
modelname = 'ConstantOnly'

if strcmp(modelname,'Covariates')
    load SurveyData
   
    
    
    NPilot = N-size(SurveyData,1);
    X = zeros(N-NPilot,size(SurveyData,2));
    X(:,1) = 1;
    X(:,2:end) = SurveyData(:,2:end);
    %X(1:NPilot,size(SurveyData,2)+1) = 1;
    XLabels = [];
    XLabels{1} = 'Constant';
    for ss = 1:(numel(SurveyLabels)-1)
        XLabels{ss+1} = SurveyLabels{1+ss};
    end
    %XLabels{(size(SurveyData,2)+1)} = 'Pilot';
end
if strcmp(modelname,'MFAC')
    load SurveyData
    NPilot = N-size(SurveyData,1);
    X = SurveyData(:,4:9);
    X = [ones(N-NPilot,1) X];
    XLabels{1} = 'Constant';
    tt = 2;
    for ss = 4:9
        XLabels{tt} = SurveyLabels{ss};
        tt = tt+1;
    end
end
if strcmp(modelname,'Sex')
    load SurveyData
    NPilot = N-size(SurveyData,1);
    X = SurveyData(:,6);
    X = [ones(N-NPilot,1) X];
    XLabels{1} = 'Constant';
    tt = 2;
    for ss = 6
        XLabels{tt} = SurveyLabels{ss};
        tt = tt+1;
    end
end
if strcmp(modelname,'ConstantOnly')~=1
    N = N-NPilot;
    % remove the pilot data
    rRD(1:NPilot) = [];
    rLT(1:NPilot) = [];
    kRD(1:NPilot) = [];
    kLT(1:NPilot) = [];
    gt(1:NPilot)  = [];
    i1(1:NPilot,:) = [];
    i2(1:NPilot,:) = [];
    ACTCOMB(1:NPilot,:) = [];
    ACTCOMBPm(1:NPilot,:,:) = [];
    D1(1:NPilot,:) = [];
    D2(1:NPilot,:) = [];
    
end


% de-mean and divide by sd for all non-binary variables
%for cc = 1:size(X,2)
%    if sum((X(:,cc)==0) + (X(:,cc)==1))~=N
%        X(:,cc) = (X(:,cc)-mean(X(:,cc)))./std(X(:,cc));
%    end
%end

XX = zeros(5*N,5*size(X,2));
XII = zeros(5*N,1);
kk = 1;
for ii = 1:N
   for vv = 1:5
       XX(5*(ii-1)+vv,((vv-1)*size(X,2)+1):(vv*size(X,2)))=X(ii,:);
       XII(kk) = ii;
       kk = kk+1;
   end
    
end


%% Initial values and setup here

nTypes = 4;

rRD = max(min(rRD,log(4)),-4); % these are used for starting values. Let's lower stupidly high ones
rLT =log(max(1-rLT,0.05));
gt = min(max(gt,-3),3);
% rLT,kLT already specified with min distance estimates

lRD = 1;
lLT = 1;

theta = zeros(N,5);
theta(:,1) = rRD;
theta(:,2) = rLT;
theta(:,3) = gt;

tmm = theta-repmat(mean(theta),[N 1]);

sigmaRRG = tmm(:,1:3)'*tmm(:,1:3)./N;
sigmaProbit = eye(2); w = sigmaProbit(2,1);
%sigma = [sigmaRRG zeros(3);zeros(3) sigmaProbit];

b = [mean(theta);zeros(size(X,2)-1,5)]';
bTrans = zeros(size(b));
bTrans(:)  = 1:numel(b);

ML = repmat(b(4:5,:)',[N 1]);

LSim = 100;

Ltype = zeros(N,2);
% latent variable for type restriction matrix
RL = zeros(nTypes-1,nTypes-1,nTypes);
for ii = 1:N
    RL = eye(2);
    if kRD(ii) == 1;
        RL(1,1) = -1;
    end
    if kLT(ii) == 1;
        RL(2,2) = -1;
    end
    Ltype(ii,:) = rmvnrnd(ML(ii,:),sigmaProbit,1,RL,[0;0]);
end
    
theta(:,4:5) = Ltype;


b = [mean(theta);zeros(size(X,2)-1,5)]';
bTrans = zeros(size(b));
bTrans(:)  = 1:numel(b);

wgrid = -0.999:0.001:0.999;

%% MH step values
MH = 30;
cMH = [0.1;0.1;0.1];
clRD = 0.1;
clLT = 0.1;

GG = 100;

%% upper and lower bounds for latent LT choices

invest = [i1 i2];
LB = invest-0.5;
UB = invest+0.5;
LB(LB<=0) = -999;
UB(UB>=100) = 999;

latentLT = zeros(N,8);

xgrid = linspace(0,100,11);

i1lb = zeros(size(i1));
i1ub = zeros(size(i1));
i1c  = zeros(size(i1));
i2lb = zeros(size(i2));
i2ub = zeros(size(i2));
i2c = zeros(size(i2));
for dd = 1:numel(i1)
    
    i1lb(dd) = xgrid(find(xgrid<=i1(dd),1,'last'));
    i1ub(dd) = xgrid(find(xgrid>=i1(dd),1,'first'));
    i2lb(dd) = xgrid(find(xgrid<=i2(dd),1,'last'));
    i2ub(dd) = xgrid(find(xgrid>=i2(dd),1,'first'));
    
end
i1c(i1ub>i1lb)=(i1(i1ub>i1lb)-i1lb(i1ub>i1lb))./(i1ub(i1ub>i1lb)-i1lb(i1ub>i1lb));
i2c(i2ub>i2lb)=(i2(i2ub>i2lb)-i2lb(i2ub>i2lb))./(i2ub(i2ub>i2lb)-i2lb(i2ub>i2lb));
%for dd = 1:numel(i1)
%    
%    i1lb(dd) = sum(xgrid<=i1lb(dd));
%    i1ub(dd) = sum(xgrid<=i1ub(dd));
%    i1lb(dd) = sum(xgrid<=i2lb(dd));
%    i2ub(dd) = sum(xgrid<=i2ub(dd));
%end
i1interp(:,:,1) = i1lb;
i1interp(:,:,2) = i1ub;
i1interp(:,:,3) = i1c;
i2interp(:,:,1) = i2lb;
i2interp(:,:,2) = i2ub;
i2interp(:,:,3) = i2c;


LTgs = 11; % grid size for interpolating lottery task likelihood
%tic
%ll = LTlogit(rLT,kLT,0.1,i1,i2,c1,c2,5,LTx0);
%toc
%ecdf(sum(ll,2))

%kill

%% priors

bPriorMU = zeros(numel(b),1);
    bPriorIV = eye(numel(b));

sigmaPriorR = eye(size(sigmaRRG));
sigmaPriorRHO = numel(size(sigmaRRG,1))+2;
sigmaProbitPriorR = eye(size(sigmaProbit));
sigmaProbitPriorRHO = numel(size(sigmaProbit,1))+2;

lRDprior = [-3,1]; % lognormal prior
lLTprior = [-4,2]; % lognormal prior


% prior distribution plots


%% Simulation starts here
trimLT = 0; % truncate i1, i2 at these values so that we don't get 1./0 problem for r>1
i1 = min(max(i1,trimLT),100-trimLT);
i2 = min(max(i2,trimLT),100-trimLT);
TypeSim = 5000;

SimSize = 20000;

thetaOUT = zeros([size(theta) SimSize]);
wOUT = zeros(SimSize,1);
sigmaRRGOUT = zeros([size(sigmaRRG) SimSize]);
bOUT = zeros([size(b) SimSize]);
typePriOUT = zeros([N 4 SimSize]);
typePrOUT  = zeros([N 4 SimSize]);
lambdaOUT = zeros([SimSize 2]);

typePr = ones(N,4)./4;

ss=1;
%%


tic
for ss = ss:SimSize
        
    % draw lRD,lLT
    l0RD = sum(RDlike(rRD,kRD,gt,lRD,ACTCOMB,ACTCOMBPm,D1,D2))-0.5*(log(lRD)-lRDprior(1)).^2./lRDprior(2);
    l0LT = sum(LTlogit(rLT,kLT,lLT,i1,i2,c1,c2,LTgs));
    for mm = 1:MH
        prop = exp(log(lRD) + clRD.*randn(1,1));
        l1 = sum(RDlike(rRD,kRD,gt,prop,ACTCOMB,ACTCOMBPm,D1,D2))-0.5*(log(prop)-lRDprior(1)).^2./lRDprior(2);
        if (l1-l0RD)>log(rand(1,1));
            l0RD=l1;
            lRD = prop;
        end
        prop = exp(log(lLT) + clLT.*randn(1,1));
        l1 = sum(LTlogit(rLT,kLT,prop,i1,i2,c1,c2,LTgs));
        if (l1-l0LT)>log(rand(1,1));
            l0LT=l1;
            lLT = prop;
        end
    end
    

    % draw elements of theta
    MU = (X*b(1:3,:)');
    for mm = 1:MH
        prop = (theta(:,1:3)' + diag(cMH)*randn(3,N))';
        tmMU = theta(:,1:3)-MU;
        pmMU = prop-MU;
        lpr0 = zeros(N,1);
        lpr1 = zeros(N,1);
        for ii = 1:N
            lpr0(ii) = tmMU(ii,:)/sigmaRRG*tmMU(ii,:)';
            lpr1(ii) = pmMU(ii,:)/sigmaRRG*pmMU(ii,:)';
        end
        l0 = RDlike(theta(:,1) ,kRD,theta(:,3),lRD,ACTCOMB,ACTCOMBPm,D1,D2) ...
            +LTlogit(theta(:,2),kLT,lLT,i1,i2,c1,c2,LTgs)...
            -0.5.*lpr0;
        l1 = RDlike(prop(:,1) ,kRD,prop(:,3),lRD,ACTCOMB,ACTCOMBPm,D1,D2) ...
            +LTlogit(prop(:,2),kLT,lLT,i1,i2,c1,c2,LTgs)...
            -0.5.*lpr1;
        
        accept = (l1-l0)>log(rand(N,1));
        theta(accept,1:3) = prop(accept,:);
        
%         %if size(X,2)==1
%             % rRD
%             [cMU,cSIGMA] = ConditionalMVN(theta(:,1:3),MU,sigmaRRG,1);
%             prop = theta(:,1)+cMH(1).*sqrt(cSIGMA).*randn(N,1);
%             l0 = RDlike(rRD ,kRD,gt,lRD,ACTCOMB,ACTCOMBPm,D1,D2) -0.5.*(cMU-rRD).^2./cSIGMA;
%             lp = RDlike(prop,kRD,gt,lRD,ACTCOMB,ACTCOMBPm,D1,D2) -0.5.*(cMU-prop).^2./cSIGMA;
%             accept = (lp-l0)>log(rand(N,1));
%             rRD(accept) = prop(accept);
%             theta(:,1) = rRD;
%             
%             % rLT
%             [cMU,cSIGMA] = ConditionalMVN(theta(:,1:3),MU,sigmaRRG,2);
%             prop = theta(:,2)+cMH(2).*sqrt(cSIGMA).*randn(N,1);            
%             l0 = LTlogit(rLT,kLT,lLT,i1,i2,c1,c2,LTgs) -0.5.*(cMU-rLT).^2./cSIGMA;
%             lp = LTlogit(prop,kLT,lLT,i1,i2,c1,c2,LTgs)-0.5.*(cMU-prop).^2./cSIGMA;
%             accept = ((lp-l0)>log(rand(N,1))); %& imag(lp1)==0;
%             rLT(accept) = prop(accept);
%             theta(:,2) = rLT;
% 
%             % gt
%             [cMU,cSIGMA] = ConditionalMVN(theta(:,1:3),MU,sigmaRRG,3);
%             prop = theta(:,3)+cMH(3).*sqrt(cSIGMA).*randn(N,1);
%             l0 = RDlike(rRD,kRD,gt,lRD,ACTCOMB,ACTCOMBPm,D1,D2)    -0.5.*(cMU-gt).^2./cSIGMA;
%             lp = RDlike(rRD,kRD,prop,lRD,ACTCOMB,ACTCOMBPm,D1,D2)  -0.5.*(cMU-prop).^2./cSIGMA;
%             accept = (lp-l0)>log(rand(N,1));
%             gt(accept) = prop(accept);
%             theta(:,3) = gt;
%             
%         %else
%         %    disp('error, need to code this bit up')
%         %    kill
%         %end
    end
    rRD = theta(:,1); rLT = theta(:,2); gt = theta(:,3);
    
    % type indicators
    % calcualte type probabilities
    TypeNorm = randn(2,TypeSim);
    ML = X*b';
    ML = ML(:,4:5);
    
    for ii = 1:N
        ts = repmat(ML(ii,:)',[1 TypeSim]) + chol(sigmaProbit,'lower')*TypeNorm;
        %tsmax = repmat(max(ts),[3 1]);
        typePr(ii,1) = mean(ts(1,:)<0 & ts(2,:)<0,2)';
        typePr(ii,2) = mean(ts(1,:)<0 & ts(2,:)>0,2)';
        typePr(ii,3) = mean(ts(1,:)>0 & ts(2,:)<0,2)';
        typePr(ii,4) = mean(ts(1,:)>0 & ts(2,:)>0,2)';
    end
    llkRD0 = RDlike(rRD,zeros(N,1),gt,lRD,ACTCOMB,ACTCOMBPm,D1,D2);
    llkRD1 = RDlike(rRD,ones(N,1),gt,lRD,ACTCOMB,ACTCOMBPm,D1,D2);
    
    llkLT0 = LTlogit(rLT,zeros(N,1),lLT,i1,i2,c1,c2,LTgs);
    llkLT1 = LTlogit(rLT,ones(N,1),lLT,i1,i2,c1,c2,LTgs);
    
    typePri(:,1) = llkRD0 + llkLT0 + log(typePr(:,1)+eps);
    typePri(:,2) = llkRD0 + llkLT1 + log(typePr(:,2)+eps);
    typePri(:,3) = llkRD1 + llkLT0 + log(typePr(:,3)+eps);
    typePri(:,4) = llkRD1 + llkLT1 + log(typePr(:,4)+eps);
    
    typePri = exp(typePri-repmat(max(typePri,[],2),[1 4]));
    typePri = typePri./repmat(sum(typePri,2),[1 4]);
    
    type = discretesampleJB(typePri);
    
    kRD = zeros(N,1); kRD(type==3 | type==4) = 1;
    kLT = zeros(N,1); kLT(type==2 | type==4) = 1;
    
    % latent variables for types
%     for ii = 1:N
%         RL = eye(2);
%         if kRD(ii) == 1;
%             RL(1,1) = -1;
%         end
%         if kLT(ii) == 1;
%             RL(2,2) = -1;
%         end
%         % avoid NaNs (not needed anymore)
%         %ok = 0;
%         %while ok ==0;
%         %    lty = rmvnrnd(ML(ii,:),sigmaProbit,1,RL,[0;0]);
%         %    ok = sum(isnan(lty))==0;
%         %end
%         
%         Ltype(ii,:) = GewekeMVTruncnormFast(ML(ii,:)',sigmaProbit,[-999;-999],[0;0],RL,1000);
%     end
%     for kk1 = 0:1
%         for kk2 = 0:1
%             % restriction matrix
%             RL = eye(2);
%             RL(1,1) = 1-2.*kk1;
%             RL(2,2) = 1-2.*kk2
%             I = (kRD == kk1) & (kLT == kk2);
%             Ltype(I,:) = GewekeMVTruncnormFast(ML(I,:)',sigmaProbit,[-999;-999],[0;0],RL,100)';
%         end
%     end
%     [kRD kLT Ltype]
%   % keep a lid on really large values
    %Ltype(abs(Ltype)==Inf) = sign(Ltype(abs(Ltype)==Inf)).*0.1.*realmax;

    Ltype(:,1) = truncnorm3(ML(:,1),1,-999.*(kRD==0),999.*(kRD==1));
    ML2cond = ML(:,2) + w.*(Ltype(:,1)-ML(:,1));
    s2cond   = (1-w.^2);
    Ltype(:,2) = truncnorm3(ML2cond,s2cond,-999.*(kLT==0),999.*(kLT==1));
    % keep a lid on really large values
    %Ltype(abs(Ltype)==Inf) = sign(Ltype(abs(Ltype)==Inf)).*0.1.*realmax;
    Ltype = min(max(Ltype,-10^6),10^-6);
    
    theta(:,4:5) = Ltype;
    
    % multinomial linear regression of theta against X
    % beta term
    sigma = zeros(5);
    sigma(1:3,1:3) = sigmaRRG;
    sigma(4:5,4:5) = sigmaProbit;
    XSinvX = zeros(numel(b));
    XSinvL = zeros(numel(b),1);
    Sinv = eye(size(b,1))/sigma;
    for ii = 1:N
        XSinvX = XSinvX + XX(XII==ii,:)'*Sinv*XX(XII==ii,:);
        XSinvL = XSinvL + XX(XII==ii,:)'*Sinv*theta(ii,:)';
    end
    D = eye(numel(b))/(XSinvX+bPriorIV);
    d = XSinvL + bPriorIV*bPriorMU;
    cholD = chol(D,'lower');
    bLong = D*d+cholD*randn(numel(b),1);
    b = bLong(bTrans);
    
    % variance terms
    E = X*b'-theta;
    
    % sigmaRRG
    ERRG = E(:,1:3);
    sigA = sigmaPriorRHO*sigmaPriorR;
    for ii = 1:N
        sigA = sigA + ERRG(ii,:)'*ERRG(ii,:);
    end
    sigInv = wishrnd(eye(size(sigA))/sigA,N+sigmaPriorRHO);
    sigmaRRG = eye(size(sigmaRRG))/sigInv;
    
    % sigmaProbit
    EProbit = E(:,4:5);
    sigA = sigmaProbitPriorRHO*sigmaProbitPriorR;
    for ii = 1:N
        sigA = sigA + EProbit(ii,:)'*EProbit(ii,:);
    end
    lpw =  (-0.5.*(sigmaProbitPriorRHO+N+2+1)).*log(1-wgrid.^2)...
           -0.5./(1-wgrid.^2).*(sigA(1,1)+sigA(2,2)-2.*wgrid.*sigA(1,2));
    pw = exp(lpw-max(lpw))./sum(exp(lpw-max(lpw)));
    w = wgrid(discretesampleJB(pw));
    sigmaProbit = [1 w;w 1];
    
    
    % store vaules
    thetaOUT(:,:,ss) = theta;
    wOUT(ss) = w;
    sigmaRRGOUT(:,:,ss) = sigmaRRG;
    bOUT(:,:,ss) = b;
    typePriOUT(:,:,ss) = typePri;
    typePrOUT(:,:,ss)  = typePr;
    lambdaOUT(ss,:) = [lRD lLT];
    
    
    disp(['Sim step ' int2str(ss) ' of ' int2str(SimSize)])
    
    disp([exp(b(1:2)') normcdf(b(3:5)');[sqrt(diag(sigmaRRG)')  w w]; [typePr(1,:) NaN]])
    
end
toc

filename = ['C:\Users\jbland\Documents\Roommates\MixModel' modelname '.mat']

eval(['save ' filename])


AnalyzeBayes

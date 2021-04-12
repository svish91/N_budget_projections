clear ;
load('NC_Bou1.mat')
%%load data
load('crIDGroups.mat');
cr=Cr_IDG;
load('InputYield2016_Apr2020_115Co_AllCrops.mat');
cate_name={'Wheat','Rice','Maize','Other Coarse Grain','Soybean','Oil Palm',...
    'Other Oil Seeds','Cotton','Sugar Crops','Fruits and Vegetable','Other Crops'};

%% Calculating NUE and Nsurplus
FAOSTAT_CoName_FAO(51) = {'Cote dIvoire'};
FAOSTAT_CoName_FAO(160) = {'Reunion'};

% fixing mongolia
%125 country ID
idx = find(ismember(FAOSTAT_CoName_115,'Mongolia')==1);
crYkgha(:,idx,32:end) = NaN;
crInkgha(:,idx,32:end) = NaN;
Nsur_allCoCate(:,idx,32:end) = NaN;
NUE_allCoCate(:,idx,32:end) = NaN;
% fixing botswana
idx = find(ismember(FAOSTAT_CoName_115,'Botswana')==1);
crYkgha(:,idx,50:end) = NaN;
crInkgha(:,idx,50:end) = NaN;
Nsur_allCoCate(:,idx,50:end) = NaN;
NUE_allCoCate(:,idx,50:end) = NaN;

NUE=crYkgha./crInkgha;
Nsur=crInkgha-crYkgha;
Ymax= NaN(11,115,6);
AICcr = NaN(11,115,6);
confInt=NaN(11,115,6,2);
YRF_Ymax = nan(11,115);
confInt_10yr = nan(11,115,2);
alpha = 0.05;

% seed
rng(100);
for idx_cr=1:11 %
for idx_co=1:115

%% fitting functions for different years
fun=@(M,x) M(1).*(x./(x+M(1)));

%%fit for first 10 years
if length(find(isnan(crInkgha(idx_cr,idx_co,1:10))==1))==10 || length(find(isnan(crYkgha(idx_cr,idx_co,1:10))==1))==10
    Ymax(idx_cr,idx_co,1) = NaN;
else
    fit10=fitnlm(reshape(crInkgha(idx_cr,idx_co,1:10),[10,1]),reshape(crYkgha(idx_cr,idx_co,1:10),[10,1]),fun,1);
    Ymax(idx_cr,idx_co,1) = fit10.Coefficients.Estimate(1);
    AICcr(idx_cr,idx_co,1) = fit10.ModelCriterion.AIC(1);
    % bootstrap
    testIn=reshape(crInkgha(idx_cr,idx_co,1:10),[10,1]);    
    testY=reshape(crYkgha(idx_cr,idx_co,1:10),[10,1]);
    confInt(idx_cr,idx_co,1,:)=f_bootstrap_YRF(testIn,testY,fun);
end
    %%fit for 20 years
if length(find(isnan(crInkgha(idx_cr,idx_co,11:20))==1))==10 || length(find(isnan(crYkgha(idx_cr,idx_co,11:20))==1))==10
    Ymax(idx_cr,idx_co,2) = NaN;
else
    fit20=fitnlm(reshape(crInkgha(idx_cr,idx_co,11:20),[10,1]),reshape(crYkgha(idx_cr,idx_co,11:20),[10,1]),fun,1);
    Ymax(idx_cr,idx_co,2) = fit20.Coefficients.Estimate(1);
    AICcr(idx_cr,idx_co,2) = fit20.ModelCriterion.AIC(1);
    % bootstrap
    testIn=reshape(crInkgha(idx_cr,idx_co,11:20),[10,1]);    
    testY=reshape(crYkgha(idx_cr,idx_co,11:20),[10,1]);
    confInt(idx_cr,idx_co,2,:)=f_bootstrap_YRF(testIn,testY,fun);

end
%%fit for 30 years
if length(find(isnan(crInkgha(idx_cr,idx_co,21:30))==1))==10 || length(find(isnan(crYkgha(idx_cr,idx_co,21:30))==1))==10
    Ymax(idx_cr,idx_co,3) = NaN;
else
    fit30=fitnlm(reshape(crInkgha(idx_cr,idx_co,21:30),[10,1]),reshape(crYkgha(idx_cr,idx_co,21:30),[10,1]),fun,1);
    Ymax(idx_cr,idx_co,3) = fit30.Coefficients.Estimate(1);
    AICcr(idx_cr,idx_co,3) = fit30.ModelCriterion.AIC(1);
    % bootstrap
    testIn=reshape(crInkgha(idx_cr,idx_co,21:30),[10,1]);    
    testY=reshape(crYkgha(idx_cr,idx_co,21:30),[10,1]);
    confInt(idx_cr,idx_co,3,:)=f_bootstrap_YRF(testIn,testY,fun);

end
%%fit for 40 years
if length(find(isnan(crInkgha(idx_cr,idx_co,31:40))==1))==10 || length(find(isnan(crYkgha(idx_cr,idx_co,31:40))==1))==10
    Ymax(idx_cr,idx_co,4) = NaN;
else
    fit40=fitnlm(reshape(crInkgha(idx_cr,idx_co,31:40),[10,1]),reshape(crYkgha(idx_cr,idx_co,31:40),[10,1]),fun,1);
    Ymax(idx_cr,idx_co,4) = fit40.Coefficients.Estimate(1);
    AICcr(idx_cr,idx_co,4) = fit40.ModelCriterion.AIC(1);
    % bootstrap
    testIn=reshape(crInkgha(idx_cr,idx_co,31:40),[10,1]);    
    testY=reshape(crYkgha(idx_cr,idx_co,31:40),[10,1]);
    confInt(idx_cr,idx_co,4,:)=f_bootstrap_YRF(testIn,testY,fun);

end
%%fit for 51 years
if length(find(isnan(crInkgha(idx_cr,idx_co,41:50))==1))==10 || length(find(isnan(crYkgha(idx_cr,idx_co,41:50))))==10
    Ymax(idx_cr,idx_co,5) = NaN;
else
    fit50=fitnlm(reshape(crInkgha(idx_cr,idx_co,41:50),[10,1]),reshape(crYkgha(idx_cr,idx_co,41:50),[10,1]),fun,1);
    Ymax(idx_cr,idx_co,5) = fit50.Coefficients.Estimate(1);
    AICcr(idx_cr,idx_co,5) = fit50.ModelCriterion.AIC(1);
    % bootstrap
    testIn=reshape(crInkgha(idx_cr,idx_co,41:50),[10,1]);    
    testY=reshape(crYkgha(idx_cr,idx_co,41:50),[10,1]);
    confInt(idx_cr,idx_co,5,:)=f_bootstrap_YRF(testIn,testY,fun);
end

%%fit for 55 years
if length(find(isnan(crInkgha(idx_cr,idx_co,51:55))==1))==5 || length(find(isnan(crYkgha(idx_cr,idx_co,51:55))==1))==5
    Ymax(idx_cr,idx_co,6) = NaN;
else
    fit55=fitnlm(reshape(crInkgha(idx_cr,idx_co,51:55),[5,1]),reshape(crYkgha(idx_cr,idx_co,51:55),[5,1]),fun,1);
    Ymax(idx_cr,idx_co,6) = fit55.Coefficients.Estimate(1);
    AICcr(idx_cr,idx_co,6) = fit55.ModelCriterion.AIC(1);
    % bootstrap
    testIn=reshape(crInkgha(idx_cr,idx_co,51:55),[5,1]);    
    testY=reshape(crYkgha(idx_cr,idx_co,51:55),[5,1]);
    confInt(idx_cr,idx_co,6,:)=f_bootstrap_YRF(testIn,testY,fun);
end
% recent 10 years
yrs = 46:55; totyrs = length(yrs);
fit55=fitnlm(reshape(crInkgha(idx_cr,idx_co,yrs),[totyrs,1]),reshape(crYkgha(idx_cr,idx_co,yrs),[totyrs,1]),fun,1);
M = fit55.Coefficients.Estimate(1);
YRF_Ymax(idx_cr,idx_co)=M;
% bootstrap
testIn=reshape(crInkgha(idx_cr,idx_co,yrs),[totyrs,1]);    
testY=reshape(crYkgha(idx_cr,idx_co,yrs),[totyrs,1]);
confInt_10yr(idx_cr,idx_co,:)=f_bootstrap_YRF(testIn,testY,fun);
end
disp([idx_cr,idx_co])
end


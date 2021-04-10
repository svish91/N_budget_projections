% program for estimating the N input for 2050 based on FAO scenario

%this program is estimating the projected N input and NUE for 2030 and 2050
%year 2006 as the baseline
clear;clc;
load('Agg_ProjectionsCrCate2050_115Co_Apr2020_AllCrops.mat')
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\iFarmData.mat','FAOSTAT_CrName_FAO')
cd('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData');
load('NC_Bou1.mat')
load('Nbudget_01-Nov-2019.mat','Nyld_kgkm_agg','Ninput_kgkm_agg','Nsur_kgkm_agg','Ninput_kgkm','Nyield_kgkm','AreaH_FAO')
Nyld_kgha_agg = Nyld_kgkm_agg./100;% km2 to ha
Ninput_kgha_agg = Ninput_kgkm_agg./100;
Nsur_kgha_agg = Nsur_kgkm_agg./100;
%%load data
load('Main_NInputYield2016_115Co_Apr2020_AllCrops.mat');
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\CropCate_AreaH_115Co_Apr2020_AllCrops.mat')

FAOSTAT_CoName_FAO(51) = {'Cote Divoire'};
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
%load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\Fert_Cr_Price_SAMData.mat')
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\CrIDGroups.mat','Cr_IDG')

cate_name={'Wheat','Rice','Maize','Other Coarse Grain','Soybean','Oil Palm',...
    'Other Oil Seeds','Cotton','Sugar Crops','Fruits and Vegetable','Other Crops'};

cd('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetWork\AGU poster project\Updating Methods 20190524\Hyperbolic Test\Uncertainty quantification\NewCoSet_115Co_Apr2020_AllCrops');
%%
clc;
co_tmp = FAOSTAT_CoName_115;

% initializing variables
ID_noProjY=nan(length(co_tmp),3);% finding Nan projected values
ID_negYmax = nan(length(co_tmp),3); % finding negative Ymax value
ID_Cases = nan(length(co_tmp),3); % finding N input projected which has been treated
ID_ubNinProj = nan(length(co_tmp),3); % countries alloted upper bound
ProjNIn2050 = nan(length(co_tmp),3);
Nsur2050 = nan(length(co_tmp),3);
NUE2050 = nan(length(co_tmp),3);
YRF_Ymax = nan(length(co_tmp),3);
bootM = nan(length(co_tmp),1000);
confInt = nan(length(co_tmp),2);
alpha = 0.05;

% NUE

Nin_t=  crInkgkm(:,:,1:55);
Ny_t = crYkgkm(:,:,1:55);
ar_t = TotAr_Cate(:,:,1:55);

numerator = reshape(nansum(Ny_t.*ar_t,1)./nansum(ar_t,1),length(co_tmp),55);
denominator = reshape(nansum(Nin_t.*ar_t,1)./nansum(ar_t,1),length(co_tmp),55);
Nue_t = reshape(numerator./denominator,length(co_tmp),55);

yrs = 51:55; totyrs = length(yrs); org_yr= 51:55;
avgNUE=nanmean(Nue_t(:,yrs),2); % estimating average NUE


% bootstrap
rng(100)
for b=1:1000
        idx=datasample(org_yr,5);
        bootM(:,b)=nanmean(Nue_t(:,idx),2);
end

for j=1:length(co_tmp)
        confInt(j,:)=prctile(bootM(j,:),[5 95]);
end
%{
%%load data
yrs = 51:55; totyrs = length(yrs); org_yr= 51:55;
avgNUE_cocr=nanmean(NUE_allCoCate(:,:,51:55),3); % estimating average NUE


yrs = 51:55; totyrs = length(yrs); %org_yr= 51:55;
avgNUE=nanmean(Nyld_kgha_agg(:,yrs)./Ninput_kgha_agg(:,yrs),2); % estimating average NUE
%}
% estimating 95th percentile of the N input
NinTmp = nanmean(Ninput_kgha_agg(:,51:55),2);
UpperLim_NIn = prctile(NinTmp,95);

%% Organizing projections to country scale
Proj_Nyield_kghacateCo2050 = nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,1)./nansum(Proj_Area_hacateCoCr2050,1);
Proj_Area_kghacateCo2050 = nansum(Proj_Area_hacateCoCr2050,1);
%%
% user input
prompt = 'Do you (A) 500 kgN or (B) upper bound as 95th percentile? ';
ubOpt = input(prompt,'s');

for idx_co = 1:length(co_tmp)
        if isnan(Proj_Nyield_kghacateCo2050(idx_co)) || Proj_Nyield_kghacateCo2050(idx_co)==0
            ProjNIn2050(idx_co,1)=0;
            Nsur2050(idx_co,2)=0;
            NUE2050(idx_co,3)=0;
             
            % finding countreis with no data
            ID_noProj(idx_co,:)= 1;
        else    
    
    %eval(['load(''C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetWork\September 2018 Work\Scenario Development\Case 5 Global Projection Scenario\AllCountries_' char(cate_name(idx_cr)) '\Co' char(FAOSTAT_CoName_FAO(idx_co)) '_' char(cate_name(idx_cr)) '.mat'')'])
        % low NUE
        ProjNIn2050(idx_co,1)= Proj_Nyield_kghacateCo2050(idx_co)/confInt(idx_co,1);
        % real data
        ProjNIn2050(idx_co,2)= Proj_Nyield_kghacateCo2050(idx_co)/avgNUE(idx_co);
        % high NUE
        ProjNIn2050(idx_co,3)= Proj_Nyield_kghacateCo2050(idx_co)/confInt(idx_co,2);
        
        if ubOpt == 'A'
            if ProjNIn2050(idx_co)>500
                ID_ubNinProj(idx_co) = 1;
                ProjNIn2050(idx_co)=500;% upper limit for N input
            end
        else
            if ProjNIn2050(idx_co,1)>UpperLim_NIn
                ID_ubNinProj(idx_co,1) = 1;
                ProjNIn2050(idx_co,1)=UpperLim_NIn;% upper limit for N input
            end
            if ProjNIn2050(idx_co,2)>UpperLim_NIn
                ID_ubNinProj(idx_co,2) = 1;
                ProjNIn2050(idx_co,2)=UpperLim_NIn;% upper limit for N input
            end
            if ProjNIn2050(idx_co,3)>UpperLim_NIn
                ID_ubNinProj(idx_co,3) = 1;
                ProjNIn2050(idx_co,3)=UpperLim_NIn;% upper limit for N input
            end
        end
        Nsur2050(idx_co,:) = ProjNIn2050(idx_co,:) - Proj_Nyield_kghacateCo2050(idx_co);
        NUE2050(idx_co,:) = Proj_Nyield_kghacateCo2050(idx_co)./ProjNIn2050(idx_co,:);
        end
end

%%
% Aggregating by country's harvested area 
Tot_ProjNIn2050Tg  = round(nansum(ProjNIn2050.*Proj_Area_kghacateCo2050'))./10^9;%round(nansum(ProjNIn2050,2).*nansum(Proj_Area_hacateCoCr2050,2))./10^9;
Tot_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCo2050.*Proj_Area_kghacateCo2050))./10^9;
Tot_NUE2050  = Tot_NYield2050Tg./Tot_ProjNIn2050Tg ;
Tot_NSur2050Tg  = Tot_ProjNIn2050Tg - Tot_NYield2050Tg;
avgNUE2050 = mean(Tot_NUE2050);




disp('Finished')
%%
% Saving the data
if ubOpt == 'A'
    eval(['save(''Results_Method1_' num2str(length(yrs)) 'yr_' num2str(500) 'ub_Feb2020_wo_CropMix_UC.mat'')']);
else
    eval(['save(''Results_Method1_' num2str(length(yrs)) 'yr_95thPub_Feb2020_wo_CropMix_115Co_UC.mat'')']);
end


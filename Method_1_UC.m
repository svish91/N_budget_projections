% program for estimating the N input for 2050 based on FAO scenario

%this program is estimating the projected N input and NUE for 2030 and 2050
%year 2006 as the baseline
clear;clc;
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\Agg_ProjectionsCrCate2050_115Co_Apr2020_AllCrops.mat')
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\iFarmData.mat','FAOSTAT_CrName_FAO')
cd('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData');
load('NC_Bou1.mat')
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\CropCate_AreaH_115Co_Apr2020_AllCrops.mat')

%%load data
load('Main_NInputYield2016_115Co_Apr2020_AllCrops.mat');

NC_Bou(170) = NaN;
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

cate_name={'Wheat','Rice','Maize','Other Coarse Grain','Soybean','Oil Palm',...
    'Other Oil Seeds','Cotton','Sugar Crops','Fruits and Vegetable','Other Crops'};

cd('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetWork\AGU poster project\Updating Methods 20190524\Hyperbolic Test\Uncertainty quantification\NewCoSet_115Co_Apr2020_AllCrops')
%%
clc;
co_tmp = FAOSTAT_CoName_115;
% initializing variables
ID_noProjY=nan(11,length(co_tmp),3);% finding Nan projected values
ID_negYmax = nan(11,length(co_tmp),3); % finding negative Ymax value
ID_Cases = nan(11,length(co_tmp),3); % finding N input projected which has been treated
ID_ubNinProj = nan(11,length(co_tmp),3); % countries alloted upper bound
ProjNIn2050 = nan(11,length(co_tmp),3);
Nsur2050 = nan(11,length(co_tmp),3);
NUE2050 = nan(11,length(co_tmp),3);
YRF_Ymax = nan(11,length(co_tmp),3);
bootM = nan(11,length(co_tmp),1000);
confInt = nan(11,length(co_tmp),2);

alpha = 0.05;
% 10 years fit
yrs = 51:55; totyrs = length(yrs); org_yr= 51:55;
avgNUE=nanmean(NUE_allCoCate(:,:,51:55),3); % estimating average NUE
% fixing Botswana
% fixing Botswana
idx = find(ismember(FAOSTAT_CoName_115,'Botswana')==1);
avgNUE(:,idx) = nanmean(NUE_allCoCate(:,idx,38:43),3);% 1998-2003 year average

% hyperbolic function to fit
%fun=@(M,x) M(1).*(x./(x+M(1)));

% estimating 95th percentile of the N input
UpperLim_NIn = NaN(11,1);
for i=1:11
    NinTmp = nanmean(reshape(crInkgha(i,:,51:55),length(co_tmp),5),2);
    UpperLim_NIn(i) = prctile(NinTmp,95);
end
%%
% bootstrap
rng(100)
for b=1:1000
        idx=datasample(org_yr,5);
        bootM(:,:,b)=nanmean(NUE_allCoCate(:,:,idx),3);
        % Botswana
        idx_B=datasample(38:43,6);
        bootM(:,25,b)=nanmean(NUE_allCoCate(:,25,idx_B),3);

end
for i=1:11
    for j=1:length(co_tmp)
        %confInt(i,j,:)=prctile(bootM(i,j,:),[100*alpha/2,100*(1-alpha/2)]);
        confInt(i,j,:)=prctile(bootM(i,j,:),[5 95]);
    end
end


%%
% user input
prompt = 'Do you (A) 500 kgN or (B) upper bound as 95th percentile? ';
ubOpt = input(prompt,'s');

for idx_cr = 1:11
    %eval(['load(''C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetWork\September 2018 Work\Scenario Development\Case 5 Global Projection Scenario\AllCountries_' char(cate_name(idx_cr)) '\Co' char(cate_name(idx_cr)) '_availData.mat'')'])
    
    for idx_co = 1:length(co_tmp)
        if isnan(Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)) || Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)==0
            ProjNIn2050(idx_cr,idx_co,:)=0;
            Nsur2050(idx_cr,idx_co,:)=0;
            NUE2050(idx_cr,idx_co,:)=0;
             
            % finding countreis with no data
            ID_noProj(idx_cr,idx_co,:)= 1;
        else    
    
    %eval(['load(''C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetWork\September 2018 Work\Scenario Development\Case 5 Global Projection Scenario\AllCountries_' char(cate_name(idx_cr)) '\Co' char(FAOSTAT_CoName_FAO(idx_co)) '_' char(cate_name(idx_cr)) '.mat'')'])
        % low NUE
        ProjNIn2050(idx_cr,idx_co,1)= Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)/confInt(idx_cr,idx_co,1);
        % real data
        ProjNIn2050(idx_cr,idx_co,2)= Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)/avgNUE(idx_cr,idx_co);
        % high NUE   
        ProjNIn2050(idx_cr,idx_co,3)= Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)/confInt(idx_cr,idx_co,2);

        if ubOpt == 'A'
            if ProjNIn2050(idx_cr,idx_co)>500
                ID_ubNinProj(idx_cr,idx_co,:) = 1;
                ProjNIn2050(idx_cr,idx_co,:)=500;% upper limit for N input
            end
        else
            if ProjNIn2050(idx_cr,idx_co,1)>UpperLim_NIn(idx_cr)
                ID_ubNinProj(idx_cr,idx_co,1) = 1;
                ProjNIn2050(idx_cr,idx_co,1)=UpperLim_NIn(idx_cr);% upper limit for N input
            end
            if ProjNIn2050(idx_cr,idx_co,2)>UpperLim_NIn(idx_cr)
                ID_ubNinProj(idx_cr,idx_co,2) = 1;
                ProjNIn2050(idx_cr,idx_co,2)=UpperLim_NIn(idx_cr);% upper limit for N input
            end
            if ProjNIn2050(idx_cr,idx_co,3)>UpperLim_NIn(idx_cr)
                ID_ubNinProj(idx_cr,idx_co,3) = 1;
                ProjNIn2050(idx_cr,idx_co,3)=UpperLim_NIn(idx_cr);% upper limit for N input
            end
        end
        Nsur2050(idx_cr,idx_co,:) = ProjNIn2050(idx_cr,idx_co,:) - Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co);
        NUE2050(idx_cr,idx_co,:) = Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)/ProjNIn2050(idx_cr,idx_co,:);
        end
    end
end % Botswana has no data in 2011-2015


% Aggregating by country's harvested area 
Tot_ProjNIn2050Tg = round(nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050),2)./10^9;
Tot_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,2))./10^9;
Tot_NUE2050  = Tot_NYield2050Tg./Tot_ProjNIn2050Tg ;
Tot_NSur2050Tg  = Tot_ProjNIn2050Tg - Tot_NYield2050Tg;
avgNUE2050 = mean(Tot_NUE2050);

% cr type
Totcr_ProjNIn2050Tg = round(nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,2),2)./10^9;
Totcr_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,2))./10^9;
Totcr_NUE2050  = Totcr_NYield2050Tg./Totcr_ProjNIn2050Tg ;
Totcr_NSur2050Tg  = Totcr_ProjNIn2050Tg - Totcr_NYield2050Tg;


%Final_estimate = table(cate_name',Tot_NYield2050Tg,Tot_ProjNIn2050Tg,Tot_NSur2050Tg,Tot_NUE2050);
% Overall sum
 % aggregate by crop type

totNin_co_kgNha = nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,1)./nansum(Proj_Area_hacateCoCr2050,1);
totNY_co_kgNha = nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,1)./nansum(Proj_Area_hacateCoCr2050,1);
totNin_co_kgNha(find(ismember(FAOSTAT_CoName_115,'Mongolia')==1)) = NaN;

totNsur_co_kgNha = totNin_co_kgNha - totNY_co_kgNha;%kg/ha
overall_totalNin_Tg = nansum(totNin_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg
overall_totalNsur_Tg = nansum(totNsur_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg
overall_totalNy_Tg = nansum(totNY_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg
%%

% overall cases 
for idx_cr =1:11
    % cases 
    PercCase_4(idx_cr) =  round(100*length(find(ID_ubNinProj(idx_cr,:) == 1))/218);
    
    % N input for these cases
    idx4 = find(ID_ubNinProj(idx_cr,:) == 1);
    NInC4(idx_cr) = nansum(ProjNIn2050(idx4).*Proj_Area_hacateCoCr2050(idx4))./10^9; %Tg
end
% percentage of cases
Overall_Cases = table(PercCase_4');
Overall_Cases.Properties.VariableNames = {'Case4'};
% N input contribution
Overall_NinputCases = table(NInC4');
Overall_NinputCases.Properties.VariableNames = {'Case4NIn'};
disp('Finished')
% Saving the data

if ubOpt == 'A'
    eval(['save(''Results_Method1_' num2str(length(yrs)) 'yr_' num2str(500) 'ub_Apr2020_115Co_UC.mat'')']);
else
    eval(['save(''Results_Method1_' num2str(length(yrs)) 'yr_95thPub_Apr2020_115Co_UC.mat'')']);
end

%% 113 countries
%{
clear;
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\Tan_Mar2020\114CountryGroup.mat')
load('Results_Method1_5yr_95thPub_Feb2020_UC.mat','ProjNIn2050','Proj_Area_hacateCoCr2050',...
    'Proj_Nyield_kghacateCoCr2050')
idx_co = Co_ID_group_X(1:113);
% Aggregating by country's harvested area 
Tot_ProjNIn2050Tg = round(nansum(ProjNIn2050(:,idx_co,:).*Proj_Area_hacateCoCr2050(:,idx_co)),2)./10^9;
Tot_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050(:,idx_co,:).*Proj_Area_hacateCoCr2050(:,idx_co),2))./10^9;
Tot_NUE2050  = Tot_NYield2050Tg./Tot_ProjNIn2050Tg ;
Tot_NSur2050Tg  = Tot_ProjNIn2050Tg - Tot_NYield2050Tg;
avgNUE2050 = mean(Tot_NUE2050);

% cr type
Totcr_ProjNIn2050Tg = round(nansum(ProjNIn2050(:,idx_co,:).*Proj_Area_hacateCoCr2050(:,idx_co),2),2)./10^9;
Totcr_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050(:,idx_co).*Proj_Area_hacateCoCr2050(:,idx_co),2))./10^9;
Totcr_NUE2050  = Totcr_NYield2050Tg./Totcr_ProjNIn2050Tg ;
Totcr_NSur2050Tg  = Totcr_ProjNIn2050Tg - Totcr_NYield2050Tg;


%Final_estimate = table(cate_name',Tot_NYield2050Tg,Tot_ProjNIn2050Tg,Tot_NSur2050Tg,Tot_NUE2050);
% Overall sum
 % aggregate by crop type

totNin_co_kgNha = nansum(ProjNIn2050(:,idx_co,:).*Proj_Area_hacateCoCr2050(:,idx_co),1)./nansum(Proj_Area_hacateCoCr2050(:,idx_co),1);
totNY_co_kgNha = nansum(Proj_Nyield_kghacateCoCr2050(:,idx_co).*Proj_Area_hacateCoCr2050(:,idx_co),1)./nansum(Proj_Area_hacateCoCr2050(:,idx_co),1);
totNin_co_kgNha(125) = NaN;

totNsur_co_kgNha = totNin_co_kgNha - totNY_co_kgNha;%kg/ha
overall_totalNin_Tg = nansum(totNin_co_kgNha.*nansum(Proj_Area_hacateCoCr2050(:,idx_co)))./10^9;%Tg
overall_totalNsur_Tg = nansum(totNsur_co_kgNha.*nansum(Proj_Area_hacateCoCr2050(:,idx_co)))./10^9;%Tg
overall_totalNy_Tg = nansum(totNY_co_kgNha.*nansum(Proj_Area_hacateCoCr2050(:,idx_co,:)))./10^9;%Tg
save('Results_Method1_5yr_95thPub_113Feb2020_UC.mat')
%}
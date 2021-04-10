clear;clc;
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\Agg_ProjectionsCrCate2050_115Co_Apr2020_AllCrops.mat')
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\Main_NInputYield2016_115Co_Apr2020_AllCrops.mat');

Ymaxs_M3 = load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetWork\AGU poster project\December 2018 Work\Ymax_all_modified_bootstrap_115Co_Apr2020_AllCrops.mat','confInt','Ymax');
Ymax = Ymaxs_M3.Ymax;
%load('Ymax_all_modified_Nov2019.mat')
Ymaxs_M2  = load('Results_Method2_10yr_95thPub_Apr2020_115Co_UC.mat','confInt','YRF_Ymax');
QNInM2 = load('Results_Method2_10yr_95thPub_Apr2020_115Co_UC.mat','ProjNIn2050','UpperLim_NIn');
YRF_Ymax = Ymaxs_M2.YRF_Ymax;
idxNa = find(YRF_Ymax ==1 );
YRF_Ymax(idxNa) = NaN;
cd('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetWork\AGU poster project\Updating Methods 20190524\Hyperbolic Test\Uncertainty quantification\NewCoSet_115Co_Apr2020_AllCrops');
%%
co_tmp = FAOSTAT_CoName_115;

% upper limit
UpperLim_NIn = QNInM2.UpperLim_NIn;
% N Input 2050 based on the fair comparison with improvement in N input for
% cases like Ymax<Yield
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\Main_NInputYield2016_115Co_Apr2020_AllCrops.mat','NUE_allCoCate');
avgNUE=nanmean(NUE_allCoCate(:,:,51:55),3);

%%%%%%%%%%%%%%%%%%%%% Estimating peojected Ymax %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear idx1 idx0 idxHigh

% Initializing variables
idx_t = [1965;1975;1985;1995;2005;2013];
screening_caseProj = NaN(11,length(co_tmp),3);
ID_Cases = NaN(11,length(co_tmp),3);
ProjNIn2050 = NaN(11,length(co_tmp),3);
Ymax2050 = NaN(11,length(co_tmp),3);
idxNegYmax = NaN(11,length(co_tmp),3);
intcept1 = NaN(11,length(co_tmp),3);
slope1 = NaN(11,length(co_tmp),3);
Rsq = NaN(11,length(co_tmp),3);
RMSE = NaN(11,length(co_tmp),3);
ID_negYmax = NaN(11,length(co_tmp),3);
ID_ubNinProj = NaN(11,length(co_tmp),3);
ID_Ymax23 = NaN(11,length(co_tmp),3);
IDquad = NaN(11,length(co_tmp),3);
Nsur2050 = NaN(11,length(co_tmp),3);
NUE2050 = NaN(11,length(co_tmp),3);
ID_noProjY = NaN(11,length(co_tmp),3);
RefYmax2050 = NaN(11,length(co_tmp),3);
% user input
prompt = 'Do you (A) 500 kgN or (B) upper bound as 95th percentile? ';
ubOpt = input(prompt,'s');
for idx_cr=1:11
    for idx_co=1:length(co_tmp)
        if isnan(Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)) || Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)==0
            ProjNIn2050(idx_cr,idx_co,:)=NaN;
            Ymax2050(idx_cr,idx_co,:)=NaN;
            Nsur2050(idx_cr,idx_co,:)=NaN;
            NUE2050(idx_cr,idx_co,:)=NaN;
            % finding countreis with no data
            ID_noProjY(idx_cr,idx_co)= 1;
        else
           % lower bound
           xx(:,1) = f_screening_Ymax(reshape(Ymaxs_M3.confInt(idx_cr,idx_co,:,1),6,1));
           % projection lower bound
           [ProjNIn2050(idx_cr,idx_co,1), screening_caseProj(idx_cr,idx_co,1), ID_ubNinProj(idx_cr,idx_co,1),Ymax2050(idx_cr,idx_co,1),...
           Nsur2050(idx_cr,idx_co,1),NUE2050(idx_cr,idx_co,1),Rsq(idx_cr,idx_co,1), RMSE(idx_cr,idx_co,1),...
           intcept1(idx_cr,idx_co,1), slope1(idx_cr,idx_co,1), ID_Cases(idx_cr,idx_co,1),ID_negYmax(idx_cr,idx_co,1),...
           RefYmax2050(idx_cr,idx_co,1),ID_Ymax23(idx_cr,idx_co,1),IDquad(idx_cr,idx_co,1)] = f_filter_projCases_lqfit(xx(:,1),QNInM2.ProjNIn2050(idx_cr,idx_co,1),...
            UpperLim_NIn(idx_cr),YRF_Ymax(idx_cr,idx_co),ubOpt, idx_t,Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), avgNUE(idx_cr,idx_co));
            % real data
            xx(:,2) = f_screening_Ymax(reshape(Ymax(idx_cr,idx_co,:),6,1));
            % projection real data
           [ProjNIn2050(idx_cr,idx_co,2), screening_caseProj(idx_cr,idx_co,2), ID_ubNinProj(idx_cr,idx_co,2),Ymax2050(idx_cr,idx_co,2),...
           Nsur2050(idx_cr,idx_co,2),NUE2050(idx_cr,idx_co,2),Rsq(idx_cr,idx_co,2), RMSE(idx_cr,idx_co,2),...
           intcept1(idx_cr,idx_co,2), slope1(idx_cr,idx_co,2), ID_Cases(idx_cr,idx_co,2),ID_negYmax(idx_cr,idx_co,2),...
           RefYmax2050(idx_cr,idx_co,2),ID_Ymax23(idx_cr,idx_co,2),IDquad(idx_cr,idx_co,2)] = f_filter_projCases_lqfit(xx(:,2),QNInM2.ProjNIn2050(idx_cr,idx_co,2),...
            UpperLim_NIn(idx_cr),YRF_Ymax(idx_cr,idx_co),ubOpt, idx_t,Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), avgNUE(idx_cr,idx_co));
            % upper bound
            xx(:,3) = f_screening_Ymax(reshape(Ymaxs_M3.confInt(idx_cr,idx_co,:,2),6,1));
           % projection lower bound
           [ProjNIn2050(idx_cr,idx_co,3), screening_caseProj(idx_cr,idx_co,3), ID_ubNinProj(idx_cr,idx_co,3),Ymax2050(idx_cr,idx_co,3),...
           Nsur2050(idx_cr,idx_co,3),NUE2050(idx_cr,idx_co,3),Rsq(idx_cr,idx_co,3), RMSE(idx_cr,idx_co,3),...
           intcept1(idx_cr,idx_co,3), slope1(idx_cr,idx_co,3), ID_Cases(idx_cr,idx_co,3),ID_negYmax(idx_cr,idx_co,3),...
           RefYmax2050(idx_cr,idx_co,3),ID_Ymax23(idx_cr,idx_co,3),IDquad(idx_cr,idx_co,3)] = f_filter_projCases_lqfit(xx(:,3),QNInM2.ProjNIn2050(idx_cr,idx_co,3),...
            UpperLim_NIn(idx_cr),YRF_Ymax(idx_cr,idx_co),ubOpt, idx_t,Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), avgNUE(idx_cr,idx_co)); 
        end
    end
end
disp('Finish')
% country
Tot_ProjNIn2050Tg  = round(nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,2))./10^9;%round(nansum(ProjNIn2050,2).*nansum(Proj_Area_hacateCoCr2050,2))./10^9;
Tot_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,2))./10^9;
Tot_NUE2050  = Tot_NYield2050Tg./Tot_ProjNIn2050Tg ;
Tot_NSur2050Tg  = Tot_ProjNIn2050Tg - Tot_NYield2050Tg;
avgNUE2050 = mean(Tot_NUE2050);

% cr type
Totcr_ProjNIn2050Tg = round(nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,2),2)./10^9;
Totcr_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,2))./10^9;
Totcr_NUE2050  = Totcr_NYield2050Tg./Totcr_ProjNIn2050Tg ;
Totcr_NSur2050Tg  = Totcr_ProjNIn2050Tg - Totcr_NYield2050Tg;

totNin_co_kgNha = nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,1)./nansum(Proj_Area_hacateCoCr2050,1);
totNY_co_kgNha = nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,1)./nansum(Proj_Area_hacateCoCr2050,1);
totNsur_co_kgNha = totNin_co_kgNha - totNY_co_kgNha;

overall_totalNin_Tg = nansum(totNin_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg
overall_totalNy_Tg = nansum(totNY_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg
overall_totalNsur_Tg = overall_totalNin_Tg - overall_totalNy_Tg;

% overall cases 
for idx_cr =1:11
    % cases 
    PercCase_1(idx_cr) = round(100*length(find(ID_Cases(idx_cr,:) == 1))/length(co_tmp));
    PercCase_2(idx_cr) = round(100*length(find(ID_Cases(idx_cr,:) == 2))/length(co_tmp));
    PercCase_3(idx_cr) =  round(100*length(find(ID_Cases(idx_cr,:) == 3))/length(co_tmp));
    PercCase_4(idx_cr) =  round(100*length(find(ID_ubNinProj(idx_cr,:) == 4))/length(co_tmp));
    PercCase_neg(idx_cr) =  round(100*length(find(ID_negYmax(idx_cr,:) == 1))/length(co_tmp));
    PercCase_noFAOYield(idx_cr) =  round(100*length(find(ID_noProjY(idx_cr,:) == 1))/length(co_tmp));
    
    % N input for these cases
    idx1 = find(ID_Cases(idx_cr,:) ==1);
    NInC1(idx_cr) = nansum(ProjNIn2050(idx1).*Proj_Area_hacateCoCr2050(idx1))./10^9; %Tg
    idx1 = find(ID_Cases(idx_cr,:) ==2);
    NInC2(idx_cr) = nansum(ProjNIn2050(idx1).*Proj_Area_hacateCoCr2050(idx1))./10^9; %Tg
    idx3 = find(ID_Cases(idx_cr,:) ==3);
    NInC3(idx_cr) = nansum(ProjNIn2050(idx3).*Proj_Area_hacateCoCr2050(idx3))./10^9; %Tg
    idx4 = find(ID_ubNinProj(idx_cr,:) ==4);
    NInC4(idx_cr) = nansum(ProjNIn2050(idx4).*Proj_Area_hacateCoCr2050(idx4))./10^9; %Tg
    idxneg = find(ID_negYmax(idx_cr,:) ==1);
    NInCneg(idx_cr) = nansum(ProjNIn2050(idxneg).*Proj_Area_hacateCoCr2050(idxneg))./10^9; %Tg
    
   
        
    % Harvested Area for these cases
    idx1 = find(ID_Cases(idx_cr,:) ==1);
    ArC1(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx1)); %ha
    idx1 = find(ID_Cases(idx_cr,:) ==2);
    ArC2(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx1)); %ha
    idx3 = find(ID_Cases(idx_cr,:) ==3);
    ArC3(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx3)); %ha
    idx4 = find(ID_ubNinProj(idx_cr,:) ==4);
    ArC4(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx4)); %ha
    idxneg = find(ID_negYmax(idx_cr,:) ==1);
    ArCneg(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idxneg)); %ha
end
% percentage of cases
Overall_Cases = table(PercCase_1',PercCase_2',PercCase_3',PercCase_4',PercCase_neg',PercCase_noFAOYield');
Overall_Cases.Properties.VariableNames = {'Case1','Case2','Case3','Case4','NegYmax','NoFAOYield2050'};
% N input contribution
Overall_NinputCases = table(NInC1',NInC2',NInC3',NInC4',NInCneg');
Overall_NinputCases.Properties.VariableNames = {'Case1NIn','Case2NIn','Case3NIn','Case4NIn','NegYmaxNIn'};
% Harvestes Area contribution
Overall_areaHCases = table(ArC1',ArC2',ArC3',ArC4',ArCneg');
Overall_areaHCases.Properties.VariableNames = {'Case1Ar','Case2Ar','Case3Ar','Case4Ar','NegYmaxAr'};

% Saving the data
if ubOpt == 'A'
    save('Results_Method3_500ub_Apr2019_lqfit_UC.mat')
else
    %save('Results_Method3_95thPub_Nov2019_lqfit_negSlope_thenAllquadCase.mat')
    save('Results_Method3_95thPub_Apr2020_lqfit_negSlope_and_UshapedQuad_UC.mat')
end

%% 113 countries
%{
clear;
load('C:\Users\svishwakarma\Documents\Research_Work\NitrogenBudgetData\Tan_Mar2020\114CountryGroup.mat')
load('Results_Method3_95thPub_Feb2020_lqfit_negSlope_and_UshapedQuad_UC.mat','ProjNIn2050','Proj_Area_hacateCoCr2050',...
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
save('Results_Method3_95thPub_113Feb2020_lqfit_negSlope_and_UshapedQuad_UC.mat')
%}
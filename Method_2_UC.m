% program for estimating the N input for 2050 based on FAO scenario

%this program is estimating the projected N input and NUE for 2030 and 2050
%year 2006 as the baseline
clear;clc;
load('Agg_ProjectionsCrCate2050_115Co_Apr2020_AllCrops.mat')
load('iFarmData.mat','FAOSTAT_CrName_FAO')
load('NC_Bou1.mat')
load('CropCate_AreaH_115Co_Apr2020_AllCrops.mat')
%%load data
load('Main_NInputYield2016_115Co_Apr2020_AllCrops.mat');
NC_Bou(170) = NaN;
co=FAOSTAT_CoName_FAO;
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

% YRF YMax
load('Ymax_10yrM2_modified_bootstrap_115Co_Apr2020_AllCrops.mat','YRF_Ymax','confInt')
load('CrIDGroups.mat','Cr_IDG')
cate_name={'Wheat','Rice','Maize','Other Coarse Grain','Soybean','Oil Palm',...
    'Other Oil Seeds','Cotton','Sugar Crops','Fruits and Vegetable','Other Crops'};

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
%YRF_Ymax = nan(11,length(co_tmp));
IDym = nan(11,length(co_tmp));

%bootM = nan(1,1000);
%bootSE = nan(1,1000);
%confInt = nan(11,length(co_tmp),2);
%alpha = 0.05;
% 10 years fit
yrs = 46:55; totyrs = length(yrs); org_yr= 51:55;
avgNUE=nanmean(NUE_allCoCate(:,:,51:55),3); % estimating average NUE

% fixing Botswana
idx = find(ismember(FAOSTAT_CoName_115,'Botswana')==1);
avgNUE(:,idx) = nanmean(NUE_allCoCate(:,idx,38:43),3);% 1998-2003 year average

% hyperbolic function to fit
fun=@(M,x) M(1).*(x./(x+M(1)));

% estimating 95th percentile of the N input
UpperLim_NIn = NaN(11,1);
for i=1:11
    NinTmp = nanmean(reshape(crInkgha(i,:,51:55),length(co_tmp),5),2);
    UpperLim_NIn(i) = prctile(NinTmp,95);
end
% user input
prompt = 'Do you (A) 500 kgN or (B) upper bound as 95th percentile? ';
ubOpt = input(prompt,'s');
%seed
rng(100)
%%
for idx_cr = 1:11
    for idx_co = 1:length(co_tmp)
        %%
        % for Mongolia and Botswana
        if idx_co == find(ismember(FAOSTAT_CoName_115,'Mongolia')==1) || ...
                (idx_co ==find(ismember(FAOSTAT_CoName_115,'Botswana')==1) && totyrs ==5)
            ProjNIn2050(idx_cr,idx_co,:)=NaN;
            Nsur2050(idx_cr,idx_co,:)=NaN;
            NUE2050(idx_cr,idx_co,:)=NaN;
        else
            % simply checking changes in Ymax
            %if M2_218.YRF_Ymax(idx_cr,FAOSTAT_CoIdx_115(idx_co)) == YRF_Ymax(idx_cr,idx_co)
           %     IDym(idx_cr,idx_co) = 1;
           % end
            % finding negative Ymax
            if YRF_Ymax(idx_cr,idx_co)<0
                ID_negYmax(idx_cr,idx_co,:) = 1;
            end
            % sometimes the projected yield is not available from FAO
            % either it is NaN or 0. Checking that condition
            if isnan(Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co)) || Proj_Nyield_kghacateCoCr2050(idx_cr,...
                    idx_co)==0 
                ProjNIn2050(idx_cr,idx_co,:)=NaN;
                Nsur2050(idx_cr,idx_co,:)=NaN;
                NUE2050(idx_cr,idx_co,:)=NaN;
                % finding countreis with no data
                ID_noProjY(idx_cr,idx_co,:)= 1;
            % M = 1 only when no data is present, so no projections will be
            % there
            elseif YRF_Ymax(idx_cr,idx_co) == 1
                ProjNIn2050(idx_cr,idx_co,:)=NaN;
                Nsur2050(idx_cr,idx_co,:)=NaN;
                NUE2050(idx_cr,idx_co,:)=NaN;
            else
                if ubOpt == 'A' % 500 kgN/ha upper limit
                    % lower bound
                    [ProjNIn2050(idx_cr,idx_co,1), Nsur2050(idx_cr,idx_co,1), NUE2050(idx_cr,idx_co,1), ID_Cases(idx_cr,idx_co,1),...
                                ID_ubNinProj(idx_cr,idx_co,1)]=f_Method2(confInt(idx_cr,idx_co,1),...
                                                          Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), 500,avgNUE(idx_cr,idx_co));
                    % real data
                    [ProjNIn2050(idx_cr,idx_co,2), Nsur2050(idx_cr,idx_co,2), NUE2050(idx_cr,idx_co,2), ID_Cases(idx_cr,idx_co,2),...
                                ID_ubNinProj(idx_cr,idx_co,2)]=f_Method2(YRF_Ymax(idx_cr,idx_co),...
                                                          Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), 500,avgNUE(idx_cr,idx_co));
                    % upper bound
                    [ProjNIn2050(idx_cr,idx_co,3), Nsur2050(idx_cr,idx_co,3), NUE2050(idx_cr,idx_co,3), ID_Cases(idx_cr,idx_co,3),...
                                ID_ubNinProj(idx_cr,idx_co,3)]=f_Method2(confInt(idx_cr,idx_co,2),...
                                                          Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), 500,avgNUE(idx_cr,idx_co));
                else % dynamic upper limit
                    % lower bound
                    [ProjNIn2050(idx_cr,idx_co,1), Nsur2050(idx_cr,idx_co,1), NUE2050(idx_cr,idx_co,1), ID_Cases(idx_cr,idx_co,1),...
                        ID_ubNinProj(idx_cr,idx_co,1)]=f_Method2(confInt(idx_cr,idx_co,1),...
                            Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), UpperLim_NIn(idx_cr),avgNUE(idx_cr,idx_co));
                    % real data
                    [ProjNIn2050(idx_cr,idx_co,2), Nsur2050(idx_cr,idx_co,2), NUE2050(idx_cr,idx_co,2), ID_Cases(idx_cr,idx_co,2),...
                        ID_ubNinProj(idx_cr,idx_co,2)]=f_Method2(YRF_Ymax(idx_cr,idx_co),...
                            Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), UpperLim_NIn(idx_cr),avgNUE(idx_cr,idx_co));
                    % upper bound
                    [ProjNIn2050(idx_cr,idx_co,3), Nsur2050(idx_cr,idx_co,3), NUE2050(idx_cr,idx_co,3), ID_Cases(idx_cr,idx_co,3),...
                        ID_ubNinProj(idx_cr,idx_co,3)]=f_Method2(confInt(idx_cr,idx_co,2),...
                            Proj_Nyield_kghacateCoCr2050(idx_cr,idx_co), UpperLim_NIn(idx_cr),avgNUE(idx_cr,idx_co));
                end
            end 
            
        end
        %%
    end
   
    
end


disp('Finished')
% Aggregating by country's harvested area 
Tot_ProjNIn2050Tg  = round(nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,2))./10^9;%round(nansum(ProjNIn2050,2).*nansum(Proj_Area_hacateCoCr2050,2))./10^9;
Tot_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,2))./10^9;
Tot_NUE2050  = Tot_NYield2050Tg./Tot_ProjNIn2050Tg ;
avgNUE2050 = mean(Tot_NUE2050);
Tot_NSur2050Tg  = Tot_ProjNIn2050Tg - Tot_NYield2050Tg;

Final_estimate = table(cate_name',Tot_NYield2050Tg,Tot_ProjNIn2050Tg,Tot_NSur2050Tg,Tot_NUE2050);


% cr type
Totcr_ProjNIn2050Tg = round(nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,2),2)./10^9;
Totcr_NYield2050Tg  = round(nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,2))./10^9;
Totcr_NUE2050  = Totcr_NYield2050Tg./Totcr_ProjNIn2050Tg ;
Totcr_NSur2050Tg  = Totcr_ProjNIn2050Tg - Totcr_NYield2050Tg;


% Overall sum
 % aggregate by crop type

totNin_co_kgNha = nansum(ProjNIn2050.*Proj_Area_hacateCoCr2050,1)./nansum(Proj_Area_hacateCoCr2050,1);
totNin_co_kgNha(find(ismember(FAOSTAT_CoName_115,'Mongolia')==1)) = NaN; % treatment for MOngolia
totNY_co_kgNha = nansum(Proj_Nyield_kghacateCoCr2050.*Proj_Area_hacateCoCr2050,1)./nansum(Proj_Area_hacateCoCr2050,1);
   

totNsur_co_kgNha = totNin_co_kgNha - totNY_co_kgNha;%kg/ha
overall_totalNin_Tg = nansum(totNin_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg
overall_totalNsur_Tg = nansum(totNsur_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg
overall_totalNy_Tg = nansum(totNY_co_kgNha.*nansum(Proj_Area_hacateCoCr2050))./10^9;%Tg

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
    idx2 = find(ID_Cases(idx_cr,:) ==2);
    NInC2(idx_cr) = nansum(ProjNIn2050(idx2).*Proj_Area_hacateCoCr2050(idx2))./10^9; %Tg
    idx3 = find(ID_Cases(idx_cr,:) ==3);
    NInC3(idx_cr) = nansum(ProjNIn2050(idx3).*Proj_Area_hacateCoCr2050(idx3))./10^9; %Tg
    idx4 = find(ID_ubNinProj(idx_cr,:) ==4);
    NInC4(idx_cr) = nansum(ProjNIn2050(idx4).*Proj_Area_hacateCoCr2050(idx4))./10^9; %Tg
    idxneg = find(ID_negYmax(idx_cr,:) ==1);
    NInCneg(idx_cr) = nansum(ProjNIn2050(idxneg).*Proj_Area_hacateCoCr2050(idxneg))./10^9; %Tg
    
    % Harvested Area for these cases
    idx1 = find(ID_Cases(idx_cr,:) ==1);
    ArC1(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx1)); %ha
    idx2 = find(ID_Cases(idx_cr,:) ==2);
    ArC2(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx2)); %ha
    idx3 = find(ID_Cases(idx_cr,:) ==3);
    ArC3(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx3)); %ha
    idx4 = find(ID_ubNinProj(idx_cr,:) ==4);
    ArC4(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idx4)); %ha
    idxneg = find(ID_negYmax(idx_cr,:) ==1);
    ArCneg(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idxneg)); %ha
    idxNoY  = find(ID_noProjY(idx_cr,:)==1);
    ArNoY(idx_cr) = nansum(Proj_Area_hacateCoCr2050(idxNoY));%ha
end
% percentage of cases
Overall_Cases = table(PercCase_1',PercCase_2',PercCase_3',PercCase_4',PercCase_neg',PercCase_noFAOYield');
Overall_Cases.Properties.VariableNames = {'Case1','Case2','Case3','Case4','NegYmax','NoFAOYield2050'};
% N input contribution
Overall_NinputCases = table(NInC1',NInC2',NInC3',NInC4',NInCneg');
Overall_NinputCases.Properties.VariableNames = {'Case1NIn','Case2NIn','Case3NIn','Case4NIn','NegYmaxNIn'};
% Harvested Area contribution
Overall_areaHCases = table(ArC1',ArC2',ArC3',ArC4',ArCneg',ArNoY');
Overall_areaHCases.Properties.VariableNames = {'Case1Ar','Case2Ar','Case3Ar','Case4Ar','NegYmaxAr','NoFAOYield2050'};

% Saving the data
if ubOpt == 'A'
    eval(['save(''Results_Method2_' num2str(length(yrs)) 'yr_' num2str(500) 'ub_Apr2020_115Co_UC.mat'')']);
else
    eval(['save(''Results_Method2_' num2str(length(yrs)) 'yr_95thPub_Apr2020_115Co_UC.mat'')']);
end

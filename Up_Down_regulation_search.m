

function [Table_positive_contra_proportional, Table_negative_contra_proportional,Table_negative_directly_proportional,Table_positive_directly_proportional] = Up_Down_regulation_search(filename_model, Objective_calc, Biomass_name)
filename_model = 'MODEL1612130000_v3_Glycerol_version_4_methanol_cultivation.xls';
m=readCbModel(filename_model);
Objective_calc="pHeme_tp";
Biomass_name="Ex_biomass";

% change of objective function to user objective 
m2=m;
objective_default=checkObjective(m2);
m2=changeObjective(m,objective_default{1},0);

% calculate max biomass value

m2=changeObjective(m2,Biomass_name,1);
max_biomass=optimizeCbModel(m2,"max").f;


% calculate min and max pHeme value
m2=changeObjective(m2,Objective_calc,1);
min_objective=optimizeCbModel(m2,"min").f;
max_objective=optimizeCbModel(m2,"max").f;

% 
m_max_all = changeRxnBounds(m2,Biomass_name, max_biomass ,'l');
m_max_all = changeRxnBounds(m2,Objective_calc, max_objective ,'l');

% Calculate Steady state for max biomass and max objective
FBA_result=optimizeCbModel(m_max_all,"max");

for i=1:length(m.rxns)
    disp (m.rxns{i})
end


% ALGORITHM to find each reaction changes 
    % Prepare m model for optimisations a) set objective to biomass

biomass_name_forced_objective="Ex_biomass";
m=changeObjective(m,objective_default{1},0);
m=changeObjective(m,biomass_name_forced_objective,1);
    % set biomass default (0;1000)
m=changeRxnBounds(m,biomass_name_forced_objective, 0 ,'l');
m=changeRxnBounds(m,biomass_name_forced_objective, 1000 ,'u');
% Check the substrate from file ...

% Set FSEOF additional constraints
Forced_objective_step = max_objective / 5 *0.9;
Forced_objective_value=0;

% Add reactions ID data
Table = table(m.rxns,'VariableNames',{'Reactions_ID'});
Table=addvars(Table,m.rxnNames,'After','Reactions_ID','NewVariableNames','Reactions_name');


Table=addvars(Table,FBA_result.x,'After','Reactions_ID');

for i=1:5
    Forced_objective_value = Forced_objective_value + Forced_objective_step;
    if Forced_objective_value<0.05
        m = changeRxnBounds(m,Objective_calc, 0.05 ,'l');
    else
        m = changeRxnBounds(m,Objective_calc, Forced_objective_value ,'l');
    end
    
    FBA_result=optimizeCbModel(m,"max");
    FBA_result_1=round(FBA_result.x,6);
    Table=addvars(Table,FBA_result_1);
        
end
Table = removevars(Table,["Var2"]);


positive_table_data=Table(:,3).Variables > 0 & Table(:,4).Variables >0 & Table(:,5).Variables >0 & Table(:,6).Variables > 0 & Table(:,6).Variables > 0;
negative_table_data=Table(:,3).Variables < 0 & Table(:,4).Variables <0 & Table(:,5).Variables <0 & Table(:,6).Variables < 0 & Table(:,6).Variables < 0;



All_numbers_table_data = Table(:,3).Variables ~=0 & Table(:,4).Variables ~=0 & Table(:,5).Variables ~=0 & Table(:,6).Variables ~= 0 & Table(:,6).Variables ~= 0;


for i=2:length(Table.Reactions_ID)
    
    Mean_table_data(i,1) = round((Table(i,3).Variables + Table(i,4).Variables + Table(i,5).Variables + Table(i,6).Variables + Table(i,7).Variables) / 5 ,6);
end

Table=addvars(Table,Mean_table_data);

for i=2:length(Table.Reactions_name)
    Weighted_mean_data_sum=0;
    Fluxes_count=0;
    for ii=1:5
            Weighted_avg_data(i,ii) = round(Table(i,ii+2).Variables / Mean_table_data(i,1) ,6);
            Weighted_mean_data_sum=Weighted_mean_data_sum + Weighted_avg_data(i,ii);
            
            if round(Table(i,ii+2).Variables,6) ~= 0
                Fluxes_count=Fluxes_count+1;
            end
    end
     Weighted_avg_data(i,6)=Weighted_mean_data_sum;
     Weighted_avg_data(i,7)=Fluxes_count;
end



Table_weighted_avg_data = array2table(Weighted_avg_data,'VariableNames',{'weighted_avg_data_1', 'weighted_avg_data_2', 'weighted_avg_data_3','weighted_avg_data_4','weighted_avg_data_5','Summ_weighted_avg_data','Count_flux_value'});

filename = 'Optimisation_result.xlsx';

% Add reaction subsystem data
Table_subsystem= cell2table(m.subSystems);
Table_subsystem.Properties.VariableNames="Subsystem";
Table=[Table,Table_subsystem];

% Add reaction stoichiometry
Table_model=readtable(filename_model,'ReadRowNames',true); 
Table=[Table,Table_model.Reaction];


% Create the final table
Table_final=addvars(Table,positive_table_data);
Table_final=addvars(Table_final,negative_table_data);
Table_final=[Table_final,Table_weighted_avg_data];


%  Perform single reaction deletion to find out essential ones.

%   Single gene deletion
[grRatio_delete, grRateKO_delete, grRateWT_delete, hasEffect_delete, delRxn_delete, fluxSolution_delete] = singleRxnDeletion(m, 'FBA');
 
 %  put the data in the table (Deletion strain growth rates (1/h))
 
  Table_essential_reaction_data = array2table(grRateKO_delete,'VariableNames',{'del_strain_GR'});
  Table_final=[Table_final,Table_essential_reaction_data];
 
 
 %  Create table with all necessary calculations
 
 writetable(Table_final,filename,'FileType','spreadsheet','sheet','all_in_one','Range','A1');

% Data Filtering 


% All dependent reactions

idx = Table_final.Mean_table_data > -100 & Table_final.Mean_table_data < 100 & Table_final.Count_flux_value >=1;
Table_dependand_reactions = Table_final(idx,:);

% write subnetwork information in to Table_final

writetable(Table_dependand_reactions,filename,'FileType','spreadsheet','sheet','Dependant_reactions','Range','A1');


% Read reaction stoichiometry info from model 


 


% All Positive correlation reactions to chosen product (HEME)
    % a) positive reactions and contra correlation
    
 Positive_contra_proportional= (Table_dependand_reactions.FBA_result_1 > Table_dependand_reactions.FBA_result_1_4 & Table_dependand_reactions.FBA_result_1_4 > 0 );
 Table_positive_contra_proportional = Table_dependand_reactions(Positive_contra_proportional,:);
 
  for i=2:length(Table_positive_contra_proportional.Reactions_ID)
    Step_weighted_precentage_data(i,1) = round((Table_positive_contra_proportional(i,1+2).Variables - Table_positive_contra_proportional(i,5+2).Variables)/Table_positive_contra_proportional(i,1+2).Variables/5*100 ,4);
 end
 
 Table_positive_contra_proportional=addvars(Table_positive_contra_proportional,Step_weighted_precentage_data);
 writetable(Table_positive_contra_proportional,filename,'FileType','spreadsheet','sheet','positive_contra_proportional','Range','A1');
 
 clear Step_weighted_precentage_data
 %   b) negative reactions and contra correlation
    
 Negative_contra_proportional= (Table_dependand_reactions.FBA_result_1 < Table_dependand_reactions.FBA_result_1_4 & Table_dependand_reactions.FBA_result_1 < 0);
 Table_negative_contra_proportional = Table_dependand_reactions(Negative_contra_proportional,:);
 
 for i=2:length(Table_negative_contra_proportional.Reactions_ID)
    Step_weighted_precentage_data(i,1) = round((Table_negative_contra_proportional(i,1+2).Variables - Table_negative_contra_proportional(i,5+2).Variables)/Table_negative_contra_proportional(i,1+2).Variables/5*100  ,4);
    
 end
 Table_negative_contra_proportional=addvars(Table_negative_contra_proportional,Step_weighted_precentage_data);
 writetable(Table_negative_contra_proportional,filename,'FileType','spreadsheet','sheet','negative_contra_proportional','Range','A1');
    
 clear Step_weighted_precentage_data
 
 
 %  c) negative reactions and directly correlation
 Negative_directly_proportional= (Table_dependand_reactions.FBA_result_1 > Table_dependand_reactions.FBA_result_1_4 & Table_dependand_reactions.FBA_result_1 < 0);
 Table_negative_directly_proportional = Table_dependand_reactions(Negative_directly_proportional,:);
 
 for i=2:length(Table_negative_directly_proportional.Reactions_ID)
    Step_weighted_precentage_data(i,1) = round((Table_negative_directly_proportional(i,5+2).Variables - Table_negative_directly_proportional(i,1+2).Variables)/ Table_negative_directly_proportional(i,5+2).Variables/5*100  ,4);
    
 end
 
 Table_negative_directly_proportional=addvars(Table_negative_directly_proportional,Step_weighted_precentage_data);
 writetable(Table_negative_directly_proportional,filename,'FileType','spreadsheet','sheet','negative_directly_proportional','Range','A1');
 
 clear Step_weighted_precentage_data
  % d) positive reactions and directly correlation
  
 Positive_directly_proportional= (Table_dependand_reactions.FBA_result_1 < Table_dependand_reactions.FBA_result_1_4 & Table_dependand_reactions.FBA_result_1 > 0);
 Table_positive_directly_proportional = Table_dependand_reactions(Positive_directly_proportional,:);
 
 for i=2:length(Table_positive_directly_proportional.Reactions_ID)
    Step_weighted_precentage_data(i,1) = round((Table_positive_directly_proportional(i,5+2).Variables - Table_positive_directly_proportional(i,1+2).Variables)/ Table_positive_directly_proportional(i,5+2).Variables/5*100  ,4);
    
 end
 
 Table_positive_directly_proportional=addvars(Table_positive_directly_proportional,Step_weighted_precentage_data);
 writetable(Table_positive_directly_proportional,filename,'FileType','spreadsheet','sheet','positive_directly_proportional','Range','A1');
end   
 
 
 
 
 
 
   

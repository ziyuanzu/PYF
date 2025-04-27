# IMPORTS
# External modules
from ppmcot import *
import copy

#model_list: the list of models calculated by the function of Get_Concretemodel_Need_Data
#target_list: the double list consist of the list of target reactions of each model, where the first reaction of each target reaction list is biomass biosynthesis reaction and the second is the reaction whose flux need to be maximized or minimized
#growth_list: the list of minimize growth rates in the function of MDF_Calculation
#growth_ratio: the list of ratios of lower bounds of growth rates to maximum growth rates
#parameter_list: the list of the total number of enzyme ratios, the upper bounds for substrate input reaction fluxes and the maximum values minus the minimum values of reaction thermodynamic driving forces
#substrate_environment: the list of the carbon source uptake reactions and the initial carbon source concentrations in environment
#metabolites_monitored: the list of metabolites delivered by the strains and the product
#metabolites_mass: the dictionary whose keys are the metabolites in metabolites_monitored and values are the molar masses of the metabolites
#concentrations_initial: the dictionary whose keys are the metabolites in metabolites_monitored and values are the initial concentrations of the metabolites in environment
#reactions_monitored: the list of dictionaries whose keys are the reactions to the uptake or release of metabolites in metabolites_monitored and values are the ratios of lower bounds of the reaction fluxes to the max reaction fluxes
#substrate_mass: the list of the molar masses of the carbon sources
#t: the number of the time points for mono-cultivation
#culture_condition: the list of the inoculum ratio of the straints and time points for co-cultivation
#map_coefficient: the list of mapping constants of the instantaneous fluxes corresponding to the mean fluxes for the mono strains
#target_metabolite: the product name
#constraints: the constraints used in the simulation, where 'kinetics and thermodynamics' consist of FBA, kinetic and thermodynamic constraints, 'kinetics' consist of FBA and kinetic constraints, 'thermodynamics' consist of FBA and thermodynamic constraints and 'FBA' consist of FBA constraint

def strategy_pyf(model_list, target_list, growth_list, growth_ratio, parameter_list, substrate_environment, metabolites_monitored, metabolites_mass, concentrations_initial, reactions_monitored, substrate_mass, t, culture_condition, target_metabolite, map_coefficient, constraints):
    number_model = len(model_list)
    monitered_production_ini = ppmcot(model_list, target_list, growth_list, growth_ratio, parameter_list, substrate_environment, metabolites_monitored, metabolites_mass, concentrations_initial, reactions_monitored, substrate_mass, t, culture_condition, map_coefficient, constraints)
    target_metabolite_production_ini = monitered_production_ini[target_metabolite]
    sensitivity = {0:{'growth_ratio':[0,0],'CS_utilization':0},1:{'growth_ratio':[0,0],'CS_utilization':0}}
    for i in range(number_model):
        growth_ratio_modified=copy.deepcopy(growth_ratio)
        map_coefficient_modified=copy.deepcopy(map_coefficient)
        growth_ratio_modified[i] -= 0.0001
        monitered_production = ppmcot(model_list, target_list, growth_list, growth_ratio_modified, parameter_list, substrate_environment, metabolites_monitored, metabolites_mass, concentrations_initial, reactions_monitored, substrate_mass, t, culture_condition, map_coefficient, constraints)
        target_metabolite_production = monitered_production[target_metabolite]
        sensitivity[i]['growth_ratio'][0]=(target_metabolite_production-target_metabolite_production_ini)/target_metabolite_production_ini/0.01
        
        growth_ratio_modified[i] += 0.0002
        monitered_production = ppmcot(model_list, target_list, growth_list, growth_ratio_modified, parameter_list, substrate_environment, metabolites_monitored, metabolites_mass, concentrations_initial, reactions_monitored, substrate_mass, t, culture_condition, map_coefficient, constraints)
        target_metabolite_production = monitered_production[target_metabolite]
        sensitivity[i]['growth_ratio'][1]=(target_metabolite_production-target_metabolite_production_ini)/target_metabolite_production_ini/0.0001
        
        
        if map_coefficient[i]<=1.99:
            map_coefficient_modified[i] += 0.01
            monitered_production = ppmcot(model_list, target_list, growth_list, growth_ratio, parameter_list, substrate_environment, metabolites_monitored, metabolites_mass, concentrations_initial, reactions_monitored, substrate_mass, t, culture_condition, map_coefficient_modified, constraints)
            target_metabolite_production = monitered_production[target_metabolite]
            sensitivity[i]['CS_utilization']=(target_metabolite_production-target_metabolite_production_ini)/target_metabolite_production_ini/0.01
    return sensitivity
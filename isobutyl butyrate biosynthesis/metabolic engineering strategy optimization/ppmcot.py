# IMPORTS
# External modules
from ETGEMs_function_ETG import *
import pandas as pd
import cobra
import ast
from cobra.io import write_sbml_model
from numpy import *
import copy
import math

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
#constraints: the constraints used in the simulation, where 'kinetics and thermodynamics' consist of FBA, kinetic and thermodynamic constraints, 'kinetics' consist of FBA and kinetic constraints, 'thermodynamics' consist of FBA and thermodynamic constraints and 'FBA' consist of FBA constraint
def ppmcot(model_list, target_list, growth_list, growth_ratio, parameter_list, substrate_environment, metabolites_monitored, metabolites_mass, concentrations_initial, reactions_monitored, substrate_mass, t, culture_condition, map_coefficient, constraints):
    deta_t = 1
    s_t = 0
    number_model = len(model_list)  
    concentrations = {}
    concentrations_change = {}
    for i in range(number_model):
        concentrations_change[i] = {}
    for mea in metabolites_monitored:
        concentrations[mea] = concentrations_initial[mea]
    
    reactions_monitored_ub = {}
    reactions_monitored_lb = {}
    for i in range(number_model):
        reactions_monitored_ub[i] = {}
        reactions_monitored_lb[i] = {}
    
    for i in range(number_model):
        for rea in reactions_monitored[i]:
            if 'reverse' in rea:
                for mea in metabolites_monitored:
                    if mea in rea:
                        reactions_monitored_ub[i][rea] = concentrations_initial[mea] * substrate_mass[i] / metabolites_mass[mea]
                        reactions_monitored_lb[i][rea] = concentrations_initial[mea] * reactions_monitored[i][rea] * substrate_mass[i] / metabolites_mass[mea]
                                
    ct = 0
    d = 0
    bn = map_coefficient[0]*culture_condition[0][0]+map_coefficient[1]*culture_condition[0][1]
    
    while ct <= culture_condition[1]-1:
        for i in range(number_model):
            model_growth = copy.deepcopy(model_list[i])
            biomass_id = target_list[i][0]
            E_total=parameter_list[i][0]
            #set the carbon source to glucose
            substrate_name=substrate_environment[i][0]
            substrate_value=parameter_list[i][1]
            biomass_value=growth_list[i]
            K_value=parameter_list[i][2]
            for rea in reactions_monitored[i]:
                if 'reverse' in rea:
                    for mea in metabolites_monitored:
                        if mea in rea:
                            reactions_monitored_ub[i][rea] = concentrations[mea] * substrate_mass[i] / metabolites_mass[mea]
                            reactions_monitored_lb[i][rea] = concentrations[mea] * reactions_monitored[i][rea] * substrate_mass[i] / metabolites_mass[mea]
                    model_growth['ub_list'][rea] = reactions_monitored_ub[i][rea] * parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i])
            if 'thermodynamics' in constraints:
                try:
                    MDF_Calculation(model_growth,biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,'gurobi')
                except:
                    pass
                else:
                    #calculate the MDF values of metabolic networks
                    B_value1=MDF_Calculation(model_growth,biomass_value,biomass_id,substrate_name,substrate_value,K_value,E_total,'gurobi')
                    #calculate the biomass yield under the MDF value
                    obj_name=target_list[i][0]
                    obj_target='maximize'
                    if constraints == 'thermodynamics':
                        max_biomass_under_mdf=Max_Growth_Rate_Calculation_thermodynamics(model_growth,obj_name,obj_target,substrate_name,substrate_value,K_value,B_value1,'gurobi')
                    elif constraints == 'kinetics and thermodynamics':
                        max_biomass_under_mdf=Max_Growth_Rate_Calculation(model_growth,obj_name,obj_target,substrate_name,substrate_value,K_value,E_total,B_value1,'gurobi')
                    model_yield = copy.deepcopy(model_list[i])
                    model_yield['lb_list'][target_list[i][0]] = max_biomass_under_mdf * growth_ratio[i]
                    for rea in reactions_monitored[i]:
                        if 'reverse' in rea:
                            for mea in metabolites_monitored:
                                if mea in rea:
                                    reactions_monitored_ub[i][rea] = concentrations[mea] * substrate_mass[i] / metabolites_mass[mea]
                                    reactions_monitored_lb[i][rea] = concentrations[mea] * reactions_monitored[i][rea] * substrate_mass[i] / metabolites_mass[mea]
                            model_yield['ub_list'][rea] = reactions_monitored_ub[i][rea] * parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i])
                            model_yield['lb_list'][rea] = reactions_monitored_lb[i][rea] * parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i])
                    E_total=parameter_list[i][0]
                    #set the carbon source to glucose
                    substrate_name=substrate_environment[i][0]
                    substrate_value1=parameter_list[i][1]
                    biomass_value1=growth_list[i]
                    K_value=parameter_list[i][2]
                    biomass_id1 = target_list[i][1]
                    try:
                        MDF_Calculation(model_yield,biomass_value1,biomass_id1,substrate_name,substrate_value1,K_value,E_total,'gurobi')
                    except:
                        pass
                    else:
                        #calculate the MDF values of metabolic networks
                        B_value2=MDF_Calculation(model_yield,biomass_value1,biomass_id1,substrate_name,substrate_value1,K_value,E_total,'gurobi')
                        if B_value2 >0: 
                            #calculate the biomass yield under the MDF value
                            obj_name1=target_list[i][1]
                            biomass_id1 = target_list[i][1]
                            B_value = B_value2
                        else:
                            obj_name1=target_list[i][0]
                            biomass_id1 = target_list[i][0]
                            B_value = B_value1
                            model_yield['ub_list'][target_list[i][1]] = 0
                        obj_target='maximize'
                        if constraints == 'thermodynamics':
                            max_biomass_under_mdf=Max_Growth_Rate_Calculation_thermodynamics(model_yield,obj_name1,obj_target,substrate_name,substrate_value1,K_value,B_value,'gurobi')
                        elif constraints == 'kinetics and thermodynamics':
                            max_biomass_under_mdf=Max_Growth_Rate_Calculation(model_yield,obj_name1,obj_target,substrate_name,substrate_value1,K_value,E_total,B_value,'gurobi')
                        biomass_value1 = max_biomass_under_mdf * 0.99
                        if constraints == 'thermodynamics':
                            [min_V,Concretemodel]=Min_Flux_Sum_Calculation_thermodynamics(model_yield,biomass_value1,biomass_id1,substrate_name,substrate_value1,K_value,B_value,'gurobi')
                        elif constraints == 'kinetics and thermodynamics':
                            [min_V,Concretemodel]=Min_Flux_Sum_Calculation(model_yield,biomass_value1,biomass_id1,substrate_name,substrate_value1,K_value,E_total,B_value,'gurobi')
                        model=model_yield['model']
                        reaction_kcat_MW=model_yield['reaction_kcat_MW']
                        reaction_g0=model_yield['reaction_g0']
                        coef_matrix=model_yield['coef_matrix']
                        metabolite_list=model_yield['metabolite_list']
                        use_result = Get_Results_Thermodynamics(model,Concretemodel,reaction_kcat_MW,reaction_g0,coef_matrix,metabolite_list)
                        for rea in reactions_monitored[i]:
                            if 'reverse' in rea:
                                for mea in metabolites_monitored:
                                    if mea in rea:
                                        concentrations_change[i][mea] = -1 * use_result['flux'][rea] * deta_t / t * bn / substrate_mass[i] * metabolites_mass[mea]/(parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i]))
                            else:
                                for mea in metabolites_monitored:
                                    if mea in rea:
                                        concentrations_change[i][mea] = use_result['flux'][rea] * deta_t / t * bn / substrate_mass[i] * metabolites_mass[mea]/(parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i]))
            else:
                obj_name=target_list[i][0]
                obj_target='maximize'
                if constraints == 'FBA':
                    max_biomass_under_mdf=Max_Growth_Rate_Calculation_FBA(model_growth,obj_name,obj_target,substrate_name,substrate_value,K_value,'gurobi')
                elif constraints == 'kinetics':
                    max_biomass_under_mdf=Max_Growth_Rate_Calculation_kinetics(model_growth,obj_name,obj_target,substrate_name,substrate_value,K_value,E_total,'gurobi')
                model_yield = copy.deepcopy(model_list[i])
                model_yield['lb_list'][target_list[i][0]] = max_biomass_under_mdf * growth_ratio[i]
                for rea in reactions_monitored[i]:
                    if 'reverse' in rea:
                        for mea in metabolites_monitored:
                            if mea in rea:
                                reactions_monitored_ub[i][rea] = concentrations[mea] * substrate_mass[i] / metabolites_mass[mea]
                                reactions_monitored_lb[i][rea] = concentrations[mea] * reactions_monitored[i][rea] * substrate_mass[i] / metabolites_mass[mea]
                        model_yield['ub_list'][rea] = reactions_monitored_ub[i][rea] * parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i])
                        model_yield['lb_list'][rea] = reactions_monitored_lb[i][rea] * parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i])
                E_total=parameter_list[i][0]
                #set the carbon source to glucose
                substrate_name=substrate_environment[i][0]
                substrate_value1=parameter_list[i][1]
                biomass_value1=growth_list[i]
                K_value=parameter_list[i][2]
                biomass_id1 = target_list[i][1]
                obj_name1=target_list[i][1]
                biomass_id1 = target_list[i][1]
                obj_target='maximize'
                if constraints == 'FBA':
                    max_biomass_under_mdf=Max_Growth_Rate_Calculation_FBA(model_yield,obj_name1,obj_target,substrate_name,substrate_value1,K_value,'gurobi')
                elif constraints == 'kinetics':
                    max_biomass_under_mdf=Max_Growth_Rate_Calculation_kinetics(model_yield,obj_name1,obj_target,substrate_name,substrate_value1,K_value,E_total,'gurobi')
                biomass_value1 = max_biomass_under_mdf * 0.99
                if constraints == 'FBA':
                    [min_V,Concretemodel]=Min_Flux_Sum_Calculation_FBA(model_yield,biomass_value1,biomass_id1,substrate_name,substrate_value1,K_value,'gurobi')
                elif constraints == 'kinetics':
                    [min_V,Concretemodel]=Min_Flux_Sum_Calculation_kinetics(model_yield,biomass_value1,biomass_id1,substrate_name,substrate_value1,K_value,E_total,'gurobi')
                model=model_yield['model']
                reaction_kcat_MW=model_yield['reaction_kcat_MW']
                reaction_g0=model_yield['reaction_g0']
                coef_matrix=model_yield['coef_matrix']
                metabolite_list=model_yield['metabolite_list']
                use_result = Get_Results_Thermodynamics(model,Concretemodel,reaction_kcat_MW,reaction_g0,coef_matrix,metabolite_list)
                for rea in reactions_monitored[i]:
                    if 'reverse' in rea:
                        for mea in metabolites_monitored:
                            if mea in rea:
                                concentrations_change[i][mea] = -1 * use_result['flux'][rea] * deta_t / t * bn / substrate_mass[i] * metabolites_mass[mea]/(parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i]))
                    else:
                        for mea in metabolites_monitored:
                            if mea in rea:
                                concentrations_change[i][mea] = use_result['flux'][rea] * deta_t / t * bn / substrate_mass[i] * metabolites_mass[mea]/(parameter_list[i][1]/(substrate_environment[i][1]*culture_condition[0][i]))
        for i in range(number_model):
            for rea in reactions_monitored[i]:
                for mea in metabolites_monitored:
                    if mea in rea:
                        concentrations[mea] += concentrations_change[i][mea]
        for mea in metabolites_monitored:
            concentrations[mea] = max(0, concentrations[mea])
        ct += deta_t
        d += 1
        bn = map_coefficient[0]*culture_condition[0][0]+map_coefficient[1]*culture_condition[0][1]
        for mea in metabolites_monitored:
            print('After {} time point, the concentration of {} is {} g/L or mM'.format(ct, mea, concentrations[mea]))
    return concentrations
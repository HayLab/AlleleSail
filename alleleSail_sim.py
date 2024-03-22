### written by Michelle Johnson and Tobin Ivy

### Import packages ###
import numpy as np
import pandas as pd
import random
import argparse
import copy
import time
import os

### Step 1: Basic Input ###
ALLELES = [['E', 'O'], ['S', 'W']]
MODS = [['germline', ['S'], 'O', ['E']]]
SEX_DET = ['autosomal', [], []]
NUM_GENS = 100
INTRO = [[1, 0, 0.2]]
D_A = [[1.0], [1.0]]
R_D = [50]
F_C = [[0, ['S'], 0.1, ['Q'], 'dominant']]
S_C = []
ADD_INTRO = []

K = 10000
N_O = 100
G_F = 10
CROSS_DICT = {}
NUM_RUNS = 10

RUN = "TEST"
FILE_NAME = "alleleSail_data/EOSW_stochastic_TEST"


### Step 2: Define Functions ###

def product(*args):
    """modified code from Python Software Foundations description of itertools' product function,
    product produces the "Cartesian product of input iterables"."""
    repeat = 1
    pools = [pool for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]

    for prod in result:
        yield prod

def product_index(*args):
    """modified code from Python Software Foundations description of itertools' product function,
    product_index tracks the indices of a given product"""
    repeat=1
    pools = [pool for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[j] for i, x in enumerate(result) for j, y in enumerate(pool)]

    for prod in result:
        yield prod


def cross(mother, father):
    """cross generates the genotypes for each offspring from the cross between mother and father"""
    
    # o_d = offspring_diploidloci, lists of all combinations of alleles for each loci
    # o_g = offspring_genotypes, all possible offspring genotypes generated from o_d
    # o_d_i = offspring_diploidloci_indices, tracks allele lineage for each locus
    # o_g_i = offspring_genotypes_indices, tracks allele lineage for each genotype
        
    # o_d = [list(product(m_a, f_a)) for f_a, m_a in zip(father, mother)]

    # o_g = list(product(*o_d))

    # o_d_i = [list(product_index(m_a, f_a)) for f_a, m_a in zip(father, mother)]

    # o_g_i = list(product(*o_d_i))
    
    o_d = [list(product(m_a, f_a)) for m_a, f_a in zip(mother, father)]

    o_g = list(product(*o_d))

    o_d_i = [list(product_index(m_a, f_a)) for m_a, f_a in zip(mother, father)]

    o_g_i = list(product(*o_d_i))
    
    
    return(o_g, o_g_i)


def allele_list(genotype, list_type):
    """allele_list generates a single list of alleles from the input genotype"""
    
    # list_type determines whether output is as a genotype or haplotypes
    # p1 determines which parent's haplotype is the primary one
    
    # g_l = genotype_list, simple list of all alleles for specified genotype
    # mo_h = mother_haplotype, list of alleles from mother
    # fa_h = father_haplotype, list of alleles from father

    if list_type == 'genotype':

        g_l = []

        for locus_index, locus in enumerate(genotype):
            for allele_index, allele in enumerate(locus):
                g_l.append(genotype[locus_index][allele_index])

        return(g_l)

    elif list_type == 'haplotype':

        mo_h = []
        fa_h = []

        for locus_index, locus in enumerate(genotype):
            mo_h.append(genotype[locus_index][0])
            fa_h.append(genotype[locus_index][1])

        return(mo_h, fa_h)


def all_option(subset, overset, replacement = 0):
    """Determines whether all elements of the subset are in the overset, either with (replacement = 1) or
    without replacement (= 0)"""
    
    subsetcopy = copy.deepcopy(subset)
    oversetcopy = copy.deepcopy(overset)
    
    if replacement:
        check = 1
        for item in subsetcopy:
            if not item in oversetcopy:
                check = 0
                break
    
        return(check)
    
    else:
        check = 1
        for item in subsetcopy:
            if item in oversetcopy:
                oversetcopy.remove(item)
                
            else:
                check = 0
                break
    
        return(check)

def fitness_cost(f_c, geno_alleles):
    """fitness_cost takes into account two distinct sources of fitness affects and combines them:
    - fitness affects from single copies of alleles that can be dominant or additive
    - fitness affects from 'genotype' affects (wherein a heterozygote might be neutral but a homozygote
    for a given allele is a lethal condition)
    
    f_c = [[sex, [allele(s)], cost, rescue, cost type], ]"""
    
    fitness = [np.zeros(len(geno_alleles[0])), np.zeros(len(geno_alleles[1]))]
    
    for cost in f_c:
        if cost[0] == 2:
            for sex in [0,1]:
                for g_i, genotype in enumerate(geno_alleles[sex]):
                    if len(cost[1]) == 1:
                        if cost[4] == 'dominant':
                            if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                                fitness[sex][g_i] += (cost[1][0] in genotype)*cost[2]
                        elif cost[4] == 'additive':
                            if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                                fitness[sex][g_i] += genotype.count(cost[1][0])*cost[2]
                            
                    elif all_option(cost[1], genotype, 0) and not any([all_option(rescue, genotype) for rescue in cost[3]]):
                        fitness[sex][g_i] += cost[2]
                        
        else:
            for g_i, genotype in enumerate(geno_alleles[cost[0]]):
                if len(cost[1]) == 1:
                    if cost[4] == 'dominant':
                        if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                            fitness[cost[0]][g_i] += (cost[1][0] in genotype)*cost[2]
                    elif cost[4] == 'additive':
                        if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                            fitness[cost[0]][g_i] += genotype.count(cost[1][0])*cost[2]

                elif all_option(cost[1], genotype, 0) and not any([all_option(rescue, genotype) for rescue in cost[3]]):
                    fitness[cost[0]][g_i] += cost[2]
                        
    for sex in [0,1]:
        for g_i, genotype in enumerate(geno_alleles[sex]):
            if fitness[sex][g_i] >= 1:
                fitness[sex][g_i] = 0
            else:
                fitness[sex][g_i] = 1 - fitness[sex][g_i]
    return(fitness)

def sterility_cost(s_c, geno_alleles):
    """sterility_cost takes into account two distinct sources of fitness affects and combines them:
    - fitness affects from single copies of alleles that can be dominant or additive
    - fitness affects from 'genotype' affects (wherein a heterozygote might be neutral but a homozygote
    for a given allele is a lethal condition)
    
    s_c = [[sex, [allele(s)], cost, rescue, cost type], ]"""
    
    fecundity = [np.zeros(len(geno_alleles[0])), np.zeros(len(geno_alleles[1]))]
    
    for cost in s_c:
        if cost[0] == 2:
            for sex in [0,1]:
                for g_i, genotype in enumerate(geno_alleles[sex]):
                    if len(cost[1]) == 1:
                        if cost[4] == 'dominant':
                            if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                                fecundity[sex][g_i] += (cost[1][0] in genotype)*cost[2]
                        elif cost[4] == 'additive':
                            if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                                fecundity[sex][g_i] += genotype.count(cost[1][0])*cost[2]
                        elif cost[4] == 'multiplicative':
                            add_sc += 'later'
                    elif all_option(cost[1], genotype, 0) and not any([all_option(rescue, genotype) for rescue in cost[3]]):
                        fecundity[sex][g_i] += cost[2]
                        
        else:
            for g_i, genotype in enumerate(geno_alleles[cost[0]]):
                if len(cost[1]) == 1:
                    if cost[4] == 'dominant':
                        if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                            fecundity[cost[0]][g_i] += (cost[1][0] in genotype)*cost[2]
                    elif cost[4] == 'additive':
                        if not any([all_option(rescue, genotype) for rescue in cost[3]]):
                            fecundity[cost[0]][g_i] += genotype.count(cost[1][0])*cost[2]
                    elif cost[4] == 'multiplicative':
                        add_sc += 'later'
                elif all_option(cost[1], genotype, 0) and not any([all_option(rescue, genotype) for rescue in cost[3]]):
                    fecundity[cost[0]][g_i] += cost[2]
                        
    for sex in [0,1]:
        for g_i, genotype in enumerate(geno_alleles[sex]):
            if fecundity[sex][g_i] >= 1:
                fecundity[sex][g_i] = 0
            else:
                fecundity[sex][g_i] = 1 - fecundity[sex][g_i]
    return(fecundity)

def drive_modification_stochastic(mod_sprin, d_b, parent_alleles, mod, mod_i, mod_dict):
    
    mod_sprin_list = [mod_sprin]
    
    for m_s_ind, mod_sprin in enumerate(mod_sprin_list): #swapped order of for loops from det model, does it matter?
        for sex in [0, 1]:
            if mod_sprin[2][mod_dict[(mod[0], mod[2], sex)]] == 0:
                mod_sprin_a_l = allele_list(mod_sprin[0], 'genotype')
                if mod[2] in mod_sprin_a_l:
                    mod_a_ind = mod_sprin_a_l.index(mod[2])
                    if ((mod_a_ind%2 == sex or 
                        (mod_a_ind%2 == sex-1 and mod_sprin_a_l[mod_a_ind+1] == mod[2])) and 
                        mod[0] == 'germline'):
                        if all_option(mod[1], parent_alleles[sex]):
                            mod_sprin_list[m_s_ind][1][mod_dict[(mod[0], mod[2], sex)]+1] = d_b[mod_i][sex][mod_dict[(mod[0], mod[2], mod[2])]]
                            mod_sprin_list[m_s_ind][2][mod_dict[(mod[0], mod[2], sex)]] = 1

                            for sub_mod_i, sub_mod in enumerate(mod[3]):
                                temp_sprin = copy.deepcopy(mod_sprin)
                                temp_sprin[0][int(mod_a_ind/2)][sex] = sub_mod
                                temp_sprin[1][mod_dict[(mod[0], mod[2], sex)]+1] = d_b[mod_i][sex][sub_mod_i]
                                mod_sprin_list.append(temp_sprin)
                                
                    elif ((mod_a_ind%2 == sex or 
                          (mod_a_ind%2 == sex-1 and mod_sprin_a_l[mod_a_ind+1] == mod[2])) and 
                          mod[0] == 'zygotic'):
                        if (all_option(mod[1][0], parent_alleles[0]) and 
                            all_option(mod[1][1], parent_alleles[1]) and 
                            all_option(mod[1][2], mod_sprin_a_l)):
                            mod_sprin_list[m_s_ind][1][mod_dict[(mod[0], mod[2], sex)]+1] = d_b[mod_i][mod_dict[(mod[0], mod[2], mod[2])]]
                            mod_sprin_list[m_s_ind][2][mod_dict[(mod[0], mod[2], sex)]] = 1

                            for sub_mod_i, sub_mod in enumerate(mod[3]):
                                temp_sprin = copy.deepcopy(mod_sprin)
                                temp_sprin[0][int(mod_a_ind/2)][sex] = sub_mod
                                temp_sprin[1][mod_dict[(mod[0], mod[2], sex)]+1] = d_b[mod_i][sub_mod_i]
                                mod_sprin_list.append(temp_sprin)
                        
                    elif ((mod_a_ind%2 == sex or 
                          (mod_a_ind%2 == sex-1 and mod_sprin_a_l[mod_a_ind+1] == mod[2])) and 
                          mod[0] == 'somatic'):
                        if all_option(mod[1], mod_sprin_a_l):
                            mod_sprin_list[m_s_ind][1][mod_dict[(mod[0], mod[2], sex)]+1] = d_b[mod_i][mod_dict[(mod[0], mod[2], mod[2])]]
                            mod_sprin_list[m_s_ind][2][mod_dict[(mod[0], mod[2], sex)]] = 1
                                                              
                            for sub_mod_i, sub_mod in enumerate(mod[3]):
                                temp_sprin = copy.deepcopy(mod_sprin)
                                temp_sprin[0][int(mod_a_ind/2)][sex] = sub_mod
                                temp_sprin[1][mod_dict[(mod[0], mod[2], sex)]+1] = d_b[mod_i][sub_mod_i]
                                mod_sprin_list.append(temp_sprin)
    
    return(mod_sprin_list)


#@profile
def stochastic_sim(alleles, mods, sex_det, 
                   num_gens, intro, d_a, r_d, f_c, s_c, add_intro, 
                   run, K, n_o, g_f, 
                   cross_dict):
    """ a function that performs a stochastic simulation, given all starting parameters
    params:
        alleles - list of locis, with each loci being a list of possible alleles at that loci. wt goes last
        mods - list of lists, each list represents a possible modification. Each modification takes the form
                [timing, [required alleles for modification], target allele, [possible replacement alleles]]
        sex_det - 
        num_gens - number of generations over which to run the simulation
        intro - introduction parameters, of the form [which sex, which genotype, introduction frequency]
        d_a - list of lists representing drive activity. each drive activity gets a list, with the values inside
                representing the frequency of success for each outcome
        r_d - list of recombination distances, values range from 0 (co-inherited) to 50 (separate chromosomes)
        f_c - list of lists representing fitness cost. Each fc takes the form 
                [sex affected, [alleles required for fc], cost, [alleles required for rescue], fc type]
        s_c - list of lists of fecundity costs (sterility cost). 
        add_intro = [[sex, genotype, release proportion, release generation, release periodicity, number of releases]], 
    returns:
        df_adults - a dataframe containing each genotype possible, and the number of adults
                of that genotype, for each generation simulated
        df_alleles - a dataframe containing each allele, and the count of that allele in the population
                for each generation simulated
        df_total - a dataframe containing total number of females and males for each generation
        cross_dict - the cross dictionary: contains each possible genotype cross, and the resulting probabilities
                of each offspring
    """
    
    # pull starting information
    num_loci = len(alleles)
    all_alleles = [allele for locus in alleles for allele in locus]
    
    diploid_loci = [list(product(allele, allele)) for allele in alleles]
    genotypes_raw = list(product(*diploid_loci))
    
    nonsense_genotypes = []
    nonsense_alleles = []
    genotypes = [[], [], [], []]
    
    if sex_det[0] == 'autosomal':
        genotypes[0] = genotypes_raw
        genotypes[1] = genotypes_raw
        for genotype in genotypes_raw:
            genotypes[2].append([a for l in genotype for a in l])
            genotypes[3].append([a for l in genotype for a in l])
    
    elif sex_det[0] == 'XY':
        for genotype in genotypes_raw:
            if "Y" in genotype[-1][0]:
                nonsense_genotypes.append(genotype)
                
            else:
                geno = [a for l in genotype for a in l]
                
                if any([all_option(sexing, geno) for sexing in sex_det[1]]):
                    genotypes[0].append(genotype)
                    genotypes[2].append(geno)
                    
                elif any([all_option(m_c, geno) for m_c in sex_det[2]]):
                    genotypes[1].append(genotype)
                    genotypes[3].append(geno)
                        
                elif 'Y' in genotype[-1][1]:
                    genotypes[1].append(genotype)
                    genotypes[3].append(geno)
                    
                else:
                    genotypes[0].append(genotype)
                    genotypes[2].append(geno)
                    
        
    elif sex_det[0] == 'ZW_viable':
        for genotype in genotypes_raw:
            geno = [a for l in genotype for a in l]
            
            # check if genotype is artifical female
            if any([all_option(sexing, geno) for sexing in sex_det[1]]):
                genotypes[0].append(genotype)
                genotypes[2].append(geno)
            # check if genotype is artifical male 
            elif any([all_option(m_c, geno) for m_c in sex_det[2]]):
                genotypes[1].append(genotype)
                genotypes[3].append(geno)
            # double W must be appended to front
            elif genotype[-1] == ['W', 'W']:
                temp_genotype = [genotype]
                temp_allele = [geno]
                genotypes[0] = temp_genotype + genotypes[0]
                genotypes[2] = temp_allele + genotypes[2]
            # W implies <2 Z's, implies natural female
            elif 'W' in genotype[-1]:
                genotypes[0].append(genotype)
                genotypes[2].append(geno)
            # no W's imply all Z's, implies natural male
            else:
                genotypes[1].append(genotype)
                genotypes[3].append(geno)

    elif sex_det[0] == 'ZW':
        for genotype in genotypes_raw:
            # we assume that WW individuals are non-viable
            if genotype[-1] == ['W', 'W']:
                nonsense_genotypes.append(genotype)
                
            else:
                geno = [a for l in genotype for a in l]

                # check if genotype is artifical female
                if any([all_option(sexing, geno) for sexing in sex_det[1]]):
                    genotypes[0].append(genotype)
                    genotypes[2].append(geno)

                # check if genotype is artifical male 
                elif any([all_option(m_c, geno) for m_c in sex_det[2]]):
                    genotypes[1].append(genotype)
                    genotypes[3].append(geno)
                
                # W implies <2 Z's, implies natural female
                elif 'W' in genotype[-1]:
                    genotypes[0].append(genotype)
                    genotypes[2].append(geno)
                
                # no W's imply all Z's, implies natural male
                else:
                    genotypes[1].append(genotype)
                    genotypes[3].append(geno)

        
    elif sex_det[0] == 'plant XY':
        # later
        later = []
        
    elif sex_det[0] == 'plant autosomal':
        # later
        later = []
        
    else:
        return('Throw error message here')
    geno_alleles = genotypes[2:]
    fitness = fitness_cost(f_c, geno_alleles)
    fertility = sterility_cost(s_c, geno_alleles)
    n_r_d = num_loci-1
    adults = [[[]], [[]]]
        
    for individual in range(int(K/2)):
        adults[0][0].append(genotypes[0][-1])
        adults[1][0].append(genotypes[1][-1])

    total_pop_list = [[],[]]
    total_pop_list[0].append(len(adults[0][0]))
    total_pop_list[1].append(len(adults[1][0]))

    for sub_intro in intro:
        for individual in range(int(sub_intro[2]*K)):
            adults[sub_intro[0]][0].append(genotypes[sub_intro[0]][sub_intro[1]])

    for a in d_a:
        a.append(1 - np.sum(a))
    d_a_copy = copy.deepcopy(d_a)
    d_b = []
    for mod_i, mod in enumerate(mods):
        if mod[0] == 'germline':
            d_b.append([d_a_copy[0], d_a_copy[1]])
            del d_a_copy[0:2]
        elif mod[0] == 'zygotic' or mod[0] == 'somatic':
            d_b.append(d_a_copy[0])
            del d_a_copy[0]
    
    # dictionary for handling cross information and how that translates to an index in the cross matrices.
    mod_dict = {}
    mod_dim_max_ind = 0

    for mod in mods:
        if (mod[0], mod[2]) in mod_dict:
            mod_dim_counter = mod_dict[(mod[0], mod[2])][0]
            
        else:
            mod_dim_counter = 0
            mod_dict[mod_dim_max_ind] = (mod[0], mod[2])
            mod_dict[(mod[0], mod[2], 0)] = mod_dim_max_ind
            mod_dim_max_ind += 1
            mod_dict[mod_dim_max_ind] = (mod[0], mod[2])
            mod_dict[(mod[0], mod[2], 1)] = mod_dim_max_ind
            mod_dim_max_ind += 1
            

        for sub_mod in mod[3]:
            mod_dict[(mod[0], mod[2], sub_mod)] = mod_dim_counter
            mod_dim_counter += 1

        # Add updated dimension counter with addition of "drive failure" sub_mod
        mod_dict[(mod[0], mod[2], mod[2])] = mod_dim_counter
        mod_dict[(mod[0], mod[2])] = mod_dim_counter + 1
    
    # Set recombination rates based on recombination distances for females and males
    if isinstance(r_d[0], int) and len(r_d) == n_r_d:
        r_d_temp = [float(r) for r in r_d]
        r_d = [r_d_temp, r_d_temp]
    elif isinstance(r_d[0], float) and len(r_d) == n_r_d:
        r_d_temp = r_d
        r_d = [r_d_temp, r_d_temp]

    elif isinstance(r_d[0], list):
        if len(r_d[0]) == len(r_d[1]) == n_r_d:
            err = 'NA'
        else:
            err = 'throw error: must have the same number of recombination rates for each sex'

    else:
        err = 'throw error: r_d must be either a single list or a pair of lists of \
        length equal to the number of recombination distances'
    
    r_d_full = []
    if n_r_d > 0:
        for r_ind, rec in enumerate(r_d):
            r_d_temp = []
            # a list comprehension for time-saving reasons
            [r_d_temp.append([0.5 + 0.5*(1-r/50), 0.5*(r/50)]) for r in rec]
            r_d_full.append(r_d_temp)

        for rec in r_d_full[0:2]:
            r_d_full.append(list(product(*rec)))

        for rec in r_d_full[2:4]:
            r_d_full.append([np.prod(rd) for rd in rec])

        r_d_full.append(list(product(*r_d_full[4:6])))
        r_d_full.append([np.prod(rd) for rd in r_d_full[6]])
    else:
        r_d_full = [[1]]*8
    
    r_d_ind_list = r_d_full[7]
    
    additonal_release_list = []
    if add_intro != []:
        for add_release in add_intro:
            additonal_release_list.append([add_release[3] + add_release[4]*g for g in range(add_release[5])])

    for gen in range(num_gens):
        adults[0].append([])
        adults[1].append([])
        
        for a_i, add_release_list in enumerate(additonal_release_list):
            if gen in add_release_list:
                for individual in range(int(add_intro[a_i][2]*K)):
                    adults[add_intro[a_i][0]][gen].append(genotypes[add_intro[a_i][0]][add_intro[a_i][1]])
                    if individual == 0:
                        print("\nintroducing", genotypes[add_intro[a_i][0]][add_intro[a_i][1]], "individual to generation", gen)
        
        if adults[1][gen] != [] and adults[0][gen] != []:
            nonsense_seen = []
            # nonsense_seen.append("oneThing")
            for mother in adults[0][gen]:
                father = random.choice(adults[1][gen])

                mother_alleles = genotypes[2][genotypes[0].index(mother)]
                father_alleles = genotypes[3][genotypes[1].index(father)]

                if (genotypes[0].index(mother), genotypes[1].index(father)) in cross_dict.keys():
                    offspring_modified = cross_dict[(genotypes[0].index(mother), genotypes[1].index(father))]

                else:
                    offspring, offspring_inds = cross(mother, father)
                    offspring_modified = [[], []]

                    for sprin_ind, sprin in enumerate(offspring):
                        sprin_inds = offspring_inds[sprin_ind]
                        r_d_ind = ''

                        if n_r_d > 0:
                            for r_event in range(n_r_d):
                                if sprin_inds[r_event][0] == sprin_inds[r_event+1][0]:
                                    r_d_ind += '0'
                                else:
                                    r_d_ind += '1'

                            for r_event in range(n_r_d):
                                if sprin_inds[r_event][1] == sprin_inds[r_event+1][1]:
                                    r_d_ind += '0'
                                else:
                                    r_d_ind += '1'

                        else:
                            r_d_ind = '0'

                        mod_sprin_list = [[sprin, 
                                           list(np.ones(mod_dim_max_ind+1)), 
                                           [0 for x in range(mod_dim_max_ind)]]]

                        mod_sprin_list[0][1][0] = r_d_ind_list[int(r_d_ind, 2)]

                        for mod_sprin in mod_sprin_list:

                            if mods != []:
                                for mod_i, mod in enumerate(mods):

                                    m_s_l = drive_modification_stochastic(mod_sprin, d_b,  
                                                                          [mother_alleles, father_alleles], 
                                                                          mod, mod_i, mod_dict)
                                    mod_sprin_list.extend(m_s_l[1:])

                        for mod_sprin in mod_sprin_list:
                            offspring_modified[0].append(mod_sprin[0])
                            offspring_modified[1].append(np.prod(mod_sprin[1]))
                    offspring_modified[1] = [x/4 for x in offspring_modified[1]]
                    cross_dict[(genotypes[0].index(mother), genotypes[1].index(father))] = offspring_modified

                n_sprin = np.random.poisson(n_o*fertility[0][genotypes[0].index(mother)
                                                            ]*fertility[1][genotypes[1].index(father)])
                new_adult_inds = np.random.choice(range(len(offspring_modified[0])), n_sprin, True, offspring_modified[1])
                new_adults = [offspring_modified[0][n_a_i] for n_a_i in new_adult_inds]

                #survival_modifier = g_f/(1+(g_f-1)*(len(adults[0][gen])+len(adults[1][gen]))/K)
                survival_modifier = g_f/(1+(g_f-1)*(total_pop_list[0][gen]+total_pop_list[1][gen])/K)

                if sex_det[0] == 'autosomal':
                    sexes = np.random.choice([0, 1], len(new_adults))
                    survival_chances = np.random.rand(len(new_adults))
                    for index, new_adult in enumerate(new_adults):
                        sex = sexes[index]
                        if survival_chances[index] <= 2/n_o*fitness[sex][genotypes[sex].index(new_adult)]*survival_modifier:
                            adults[sex][gen+1].append(new_adult)
                        # if np.random.binomial(1, 2/n_o*fitness[sex][genotypes[sex].index(new_adult)]*survival_modifier):
                        #     adults[sex][gen+1].append(new_adult)

                else: # (sex_det[0] == 'XY') or (sex_det[0] == 'ZW'): # or ZW_viable!!
                    survival_chances = np.random.rand(len(new_adults))
                    for survival_index, new_adult in enumerate(new_adults):
                        # new_adult_a = allele_list(new_adult, 'genotype')
                        #if any([all_option(s_c, new_adult_a) for s_c in sex_det[2]]):
                        if new_adult in genotypes[1]:
                            #if np.random.binomial(1, 2/n_o*fitness[1][genotypes[1].index(new_adult)]*survival_modifier):
                            if survival_chances[survival_index] <= 2/n_o*fitness[1][genotypes[1].index(new_adult)]*survival_modifier:
                                adults[1][gen+1].append(new_adult)
                        elif new_adult in genotypes[0]: #any([all_option(s_c, new_adult_a) for s_c in sex_det[1]]):
                            #if np.random.binomial(1, 2/n_o*fitness[0][genotypes[0].index(new_adult)]*survival_modifier):
                            if survival_chances[survival_index] <= 2/n_o*fitness[0][genotypes[0].index(new_adult)]*survival_modifier:
                                adults[0][gen+1].append(new_adult)
                        else:
                            err = 'throw error here'
                            nonsense_seen.append(str(new_adult))
                            nonsense_seen.append("a")

            total_pop_list[0].append(len(adults[0][gen+1]))
            total_pop_list[1].append(len(adults[1][gen+1]))
    
    if nonsense_seen != []:
        print("Nonsense Genotype alert for: ")
        print(set(nonsense_seen))
    
    adults_temp = [[], []]

    for sex in [0, 1]:
        for n in range(num_gens+1):
            gen_count = []
            ClvR_count = []
            non_ClvR_count = []
            for geno in genotypes[sex]:
                gen_count.append(adults[sex][n].count(geno))
            adults_temp[sex].append(gen_count)
        adults_temp[sex] = np.array([np.array(x) for x in adults_temp[sex]])
    
    # generate df of total population counts
    df_total_females = pd.DataFrame(np.sum(adults_temp[0], axis = 1), columns = ['Count'])
    df_total_females['Generation'] = range(num_gens+1)
    df_total_females['Sex'] = 'Female'
    
    df_total_males = pd.DataFrame(np.sum(adults_temp[1], axis = 1), columns = ['Count'])
    df_total_males['Generation'] = range(num_gens+1)
    df_total_males['Sex'] = 'Male'
    df_total = pd.concat([df_total_females, df_total_males])
    df_total = df_total.reset_index(drop=True)
    df_total['Run'] = run

    #generate df of total population - NOT including those that are introduced
    tpc_final = [[],[]]
    len_tpc = len(total_pop_list[0])
    for n in range(num_gens+1):
        if len_tpc > n:
            tpc_final[0].append(total_pop_list[0][n])
            tpc_final[1].append(total_pop_list[1][n])
        else:
            tpc_final[0].append(0)
            tpc_final[1].append(0)
            
    df_total_pop_f = pd.DataFrame(tpc_final[0], columns = ['Count'])
    df_total_pop_f['Generation'] = range(num_gens+1)
    df_total_pop_f['Sex'] = 'Female'

    df_total_pop_m = pd.DataFrame(tpc_final[1], columns = ['Count'])
    df_total_pop_m['Generation'] = range(num_gens+1)
    df_total_pop_m['Sex'] = 'Male'
    df_total_pop = pd.concat([df_total_pop_f, df_total_pop_m])
    df_total_pop = df_total_pop.reset_index(drop=True)
    df_total_pop['Run'] = run
    
    # generate df of genotypes (called 'adults')
    df_females = pd.DataFrame(adults_temp[0], columns = [str(geno) for geno in genotypes[0]])
    df_females['Generation'] = range(num_gens+1)
    df_females = df_females.melt(id_vars = 'Generation', var_name = 'Genotype', value_name = 'Count')
    df_females['Sex'] = 'Female'

    df_males = pd.DataFrame(adults_temp[1], columns = [str(geno) for geno in genotypes[1]])
    df_males['Generation'] = range(num_gens+1)
    df_males = df_males.melt(id_vars = 'Generation', var_name = 'Genotype', value_name = 'Count')
    df_males['Sex'] = 'Male'
    df_adults = pd.concat([df_females, df_males])
    df_adults = df_adults.reset_index(drop=True)
    df_adults['Run'] = run

    # generate genotype-bearing allele count
    alleles_temp = [[], []]
    for sex in [0, 1]:
        for allele in all_alleles:
            a_count = np.zeros(num_gens+1)
            for geno_ind, geno in enumerate(genotypes[sex+2]):
                if allele in geno:
                    a_count += adults_temp[sex][:,geno_ind]
            alleles_temp[sex].append(a_count)
        alleles_temp[sex] = np.array([np.array(x) for x in alleles_temp[sex]])

    df_females_a = pd.DataFrame(np.transpose(alleles_temp[0]), columns = [allele for allele in all_alleles])
    df_females_a['Generation'] = range(num_gens+1)
    df_females_a = df_females_a.melt(id_vars = 'Generation', var_name = 'Allele', value_name = 'Count')
    df_females_a['Sex'] = 'Female'

    df_males_a = pd.DataFrame(np.transpose(alleles_temp[1]), columns = [allele for allele in all_alleles])
    df_males_a['Generation'] = range(num_gens+1)
    df_males_a = df_males_a.melt(id_vars = 'Generation', var_name = 'Allele', value_name = 'Count')
    df_males_a['Sex'] = 'Male'
    df_alleles = pd.concat([df_females_a, df_males_a])
    df_alleles = df_alleles.reset_index(drop=True)
    df_alleles.drop(df_alleles.index[(df_alleles["Sex"] == 'Female') & (df_alleles["Allele"] == 'Y')].tolist(), 
                    inplace = True)
    df_alleles['Run'] = run

    # generate df of PURE allele counts
    new_alleles_temp = [[], []]
    for sex in [0, 1]:
        for allele in all_alleles:
            a_count = [0]*(num_gens+1)
            for geno_index, geno in enumerate(genotypes[sex+2]):
                a_count += adults_temp[sex][:,geno_index] * geno.count(allele)
            new_alleles_temp[sex].append(a_count)
        new_alleles_temp[sex] = np.array([np.array(x) for x in new_alleles_temp[sex]])

    new_df_females_a = pd.DataFrame(np.transpose(new_alleles_temp[0]), columns = [allele for allele in all_alleles])
    new_df_females_a['Generation'] = range(num_gens+1)
    new_df_females_a = new_df_females_a.melt(id_vars = 'Generation', var_name = 'Allele', value_name = 'Count')
    new_df_females_a['Sex'] = 'Female'

    new_df_males_a = pd.DataFrame(np.transpose(new_alleles_temp[1]), columns = [allele for allele in all_alleles])
    new_df_males_a['Generation'] = range(num_gens+1)
    new_df_males_a = new_df_males_a.melt(id_vars = 'Generation', var_name = 'Allele', value_name = 'Count')
    new_df_males_a['Sex'] = 'Male'
    new_df_alleles = pd.concat([new_df_females_a, new_df_males_a])
    new_df_alleles = new_df_alleles.reset_index(drop=True)
    new_df_alleles.drop(new_df_alleles.index[(new_df_alleles["Sex"] == 'Female') & (new_df_alleles["Allele"] == 'Y')].tolist(), 
                    inplace = True)
    new_df_alleles['Run'] = run
    
    return(df_adults, df_alleles, df_total, new_df_alleles, df_total_pop, cross_dict)


def run_stochastic_sim(alleles = ALLELES, mods = MODS, sex_det = SEX_DET,
                       num_gens = NUM_GENS, intro = INTRO, d_a = D_A, r_d = R_D,
                       f_c = F_C, s_c = S_C, add_intro = ADD_INTRO, run = RUN,
                       k = K, n_o = N_O, g_f = G_F, cross_dict_orig = CROSS_DICT,
                       num_runs = NUM_RUNS, file_name = FILE_NAME):
    
    cross_dict = copy.deepcopy(cross_dict_orig)

    for i in range(1, num_runs+1):
        run_label = run + "_" + str(i)
        startTime = time.time()
        simulation = stochastic_sim(alleles, mods, sex_det, 
                              num_gens, intro, d_a, r_d, f_c, s_c, add_intro, 
                              run_label, k, n_o, g_f, 
                              cross_dict)
        endTime = time.time()
        
        for file_i, type in enumerate(["_genotype", "_allele", "_total", "_NEWallele", "_total_pop"]):
            df = simulation[file_i]
            fn = file_name + str(type) + ".csv"

            df.to_csv(fn, sep=',', mode = 'a')

        cross_dict = simulation[5]
        
        print("time taken on run " + str(i) + ": " + str(endTime - startTime))

    print(str(run) + " appended to " + str(file_name) + "\n")

def main():
    global ALLELES, MODS, SEX_DET, NUM_GENS, INTRO, D_A, R_D, F_C, S_C, ADD_INTRO
    global K, N_O, G_F, CROSS_DICT, NUM_RUNS, RUN, FILE_NAME

    # initialize argument parser
    #formatter = lambda prog: argparse.HelpFormatter(prog, width=os.get_terminal_size()[0])
    parser = argparse.ArgumentParser(#formatter_class=formatter,
                                     description="run a stochastic simulation of allele sail behavior")

    # add all possible arguments
    parser.add_argument("-a", "--alleles", help = "alleles being used for simulation",
                        nargs = '+', default = ALLELES, type=list)
    parser.add_argument("-m", "--mods", help = "modifications being used. Honestly? Never touch this",
                        type=list, default=MODS)
    parser.add_argument("--somatic", help = "sail cleavage occurs somaticly, as opposed to in the germline",
                        action="store_true")
    
    parser.add_argument("-sex", "--sex_determination", 
                        help = "which sex determination system to use, for suppression only. \
                        Options include XY, XYsterile, ZW, ZWsterile, ZW_viable, and ZW_Z_Wind",
                        type=str, default = "modification")

    # parser.add_argument("--XY", help = "use an XY determination, where our sail (S) cuts aromatase",
    #                     action="store_true")
    # parser.add_argument("--XYsterile", help = "use an XY determination, where our sail (S) cuts aromatase AND a fertility gene",
    #                     action="store_true")
    # parser.add_argument("--ZW", help = "use a ZW determination, where our sail (S) cuts aromatase",
    #                     action="store_true")
    # parser.add_argument("--ZWsterile", help = "use an ZW determination, where our sail (S) cuts aromatase AND a fertility gene",
    #                     action="store_true")
    # parser.add_argument("--ZW_viable", help = "use a ZW determination, where WW is viable and our sail (S) cuts aromatase",
    #                     action="store_true")
    
    parser.add_argument("--MC", help = "also have maternal carryover, where our sail (S) cuts aromatase",
                        action="store_true")
    parser.add_argument("-n", "--num_gens", help = "number of generations to run simulation for. \
                        Default is 100", type=int, default=NUM_GENS)
    parser.add_argument("-i", "--intro", help = "introduction frequency of males homozygous",
                        type = float, default=0.2)
    parser.add_argument("-e", "--efficiency", help = "float cleavage efficiency of our sail, in males and females. \
                        Default is 100%% (1.0) cleavage in males and females. Takes the format -e female -e male", 
                        nargs = 1, type=float, action = "append", default=[])
    parser.add_argument("-rd", "--recomb_dist", help = "list of recombination distances between loci. \
                        Defaule it R_D = [50], for maximum recomb. distance between loci 1 and 2. For more loci, \
                        must include more recomb_distances", type = int,
                        nargs = '+', default=R_D)
    parser.add_argument("-fc", "--fitness_cost", help = "list of fitness costs. structure is sex, \
                        alleles that incur the cost, the fitness cost itself. ex: -fc 2 EE 0.2 \
                        Default is no costs", nargs = '+', action = 'append',
                        default=[])
    parser.add_argument("--dominant_FC", help = "make fitness costs dominant, as opposed to additive",
                        action="store_true")
    parser.add_argument("-sc", "--sterility_cost", help = "list of sterility costs. \
                        Default is no sc, takes the same format as -fc", nargs = '+', action = 'append',
                        default=[])
    parser.add_argument("-ai", "--additional_intros", help = "list of additional releases. Takes the format" \
                        " [[sex, genotype, intro_percent, first release, timing between releases,"  \
                        " number of total additional releases]]", nargs = '+', action = 'append',
                        default=[])
    # note: run_label_type is required
    parser.add_argument("-l", "--run_label_type", help = "type of label for the runs, generally a combination of" \
                        " fitness costs and intro frequencies, for labeling purposes. Will automatically add" \
                        " introduction frequency and fitness cost, for inputs 'intro_frequency' and 'fitCosts'",
                        default=RUN, nargs = '+') #, choices=["intro_frequency", "fitCosts"])
    parser.add_argument("-k", help = "integer carrying capacity", type = int, default=K)
    parser.add_argument("-no", "--num_offspring", help = "integer number of offspring, \
                        Default 100", type = int, default=N_O)
    parser.add_argument("-gf", "--growth_factor", help="integer growth factor, default is 10",
                        type = int, default=G_F)
    parser.add_argument("-r", "--num_runs", help = "integer number of runs, default is 10",
                        type = int, default=NUM_RUNS)
    parser.add_argument("-fn", "--file_name", help="file name for data to be appended to", default=FILE_NAME)

    args = parser.parse_args()

    print(args.recomb_dist)
    exit

    if not args.efficiency:
        args.efficiency = D_A

    # handle type of simulation
    sex_det = ['autosomal', [], []]

    if args.sex_determination == "XY":
        args.alleles = [['E', 'O'], ['S', 'W'], ['X', 'Y']]
        sex_det = ['XY', [], [['E', 'E']]]
        if args.recomb_dist == R_D:
            args.recomb_dist = [50, 50]

    elif args.sex_determination == "XYsterile":
        args.alleles = [['I', 'F'], ['E', 'O'], ['S', 'W'], ['X', 'Y']]
        sex_det = ['XY', [], [['E', 'E']]]
        args.mods = [['germline', ['S'], 'O', ['E']], ['germline', ['S'], 'F', ['I']]]
        args.sterility_cost.append([0, "II", 1.0])
        if args.recomb_dist == R_D:
            args.recomb_dist = [50, 50]

    elif args.sex_determination == "ZW":
        args.alleles = [['E', 'O'], ['S', 'D'], ['Z', 'W']]
        sex_det = ['ZW', [], [['E', 'E']]]
        if args.recomb_dist == R_D:
            args.recomb_dist = [50, 50]

    elif args.sex_determination == "ZWsterile":
        args.alleles = [['I', 'F'], ['E', 'O'], ['S', 'D'], ['Z', 'W']]
        sex_det = ['ZW', [], [['E', 'E']]]
        args.mods = [['germline', ['S'], 'O', ['E']], ['germline', ['S'], 'F', ['I']]]
        args.sterility_cost.append([0, "II", 1.0])
        if args.recomb_dist == R_D:
            args.recomb_dist = [50, 50]

    elif args.sex_determination == "ZW_viable":
        args.alleles = [['E', 'O'], ['S', 'D'], ['Z', 'W']]
        sex_det = ['ZW_viable', [], [['E', 'E']]]
        if args.recomb_dist == R_D:
            args.recomb_dist = [50, 50]

    elif args.sex_determination == "ZW_Z_Wind":
        args.alleles = [['E', 'O'], ['Z*', 'Z', 'W']]
        sex_det = ['ZW', [], [['E', 'E']]]
        args.mods = [['germline', ['Z*'], 'O', ['E']]]

    if args.somatic:
        new_mod_list = []
        for mod in args.mods:
            new_mod_list.append(['somatic', mod[1], mod[2], mod[3]])
        args.mods = new_mod_list

    if args.MC:
        new_arg = ['zygotic', [['S'], [], []], 'O', ['E']]
        if args.mods[0][0] == 'somatic':
            args.mods = [new_arg, args.mods[0]]
        else:
            args.mods.append(new_arg)
            
        # also, now that there are two mods there are two efficiencies
        if args.efficiency == D_A:
            args.efficiency = [[1.0], [1.0], [1.0]] # female germline, male germline, zygotic
   
        

    run_label = []

    ## reformat introduction
    args.intro = [[1, 0, args.intro]]

    fitness_type = "additive"
    if args.dominant_FC:
        fitness_type = "dominant"

    ## handle fitness costs
    fitness_cost = []
    for sex, allele, cost in args.fitness_cost:
        fitness_cost.append([int(sex), list(allele), float(cost), ['Q'], fitness_type])

    ## handle fitness costs
    fitness_cost = []
    for sex, allele, cost in args.fitness_cost:
        fitness_cost.append([int(sex), list(allele), float(cost), ['Q'], fitness_type])
    
    print(fitness_cost)

    sterility_cost = []
    for sex, allele, cost in args.sterility_cost:
        sterility_cost.append([int(sex), list(allele), float(cost), ['Q'], "additive"])
    
    print(sterility_cost)

    ## handle additional intros
    add_intros = []
    print(args.additional_intros)
    for sex, genotype, intro_percent, release_1, gap, repeats in args.additional_intros:
        add_intros.append([int(sex), int(genotype), float(intro_percent), 
                           int(release_1), int(gap), int(repeats)])

    # if we have inputs for run_label_type, we must add all of them
    if type(args.run_label_type) == list:
        # check option 1
        if "intro_frequency" in args.run_label_type:
            run_label.append("introFreq_" + str(args.intro[0][2]))
            args.run_label_type.remove("intro_frequency")
        # check option 2
        if "fitCosts" in args.run_label_type:
            if args.fitness_cost == []:
                run_label.append("cost_0")
            for cost_list in args.fitness_cost:
            # [[0, ['S'], 0.1, ['Q'], 'dominant']]
                run_label.append(cost_list[1][0] + "cost_" + str(cost_list[2]))
            args.run_label_type.remove("fitCosts")
        # add additional labels
        run_label.extend(args.run_label_type)
        # combine all options
        run_label = "_".join(run_label)
    else:
        run_label = RUN

        
    cross_dict = {}

    ## unpack vars
    alleles, mods, somaticFlag, sex_determination, MCflag, num_gens, \
        intro, efficiency, recomb_dist, old_fitness_cost, dominantFCFlag, old_sterility_cost, old_add_intros, \
        run_label_type, k, num_offspring, \
        growth_factor, num_runs, file_name = vars(args).values()

    ## run simulation !
    run_stochastic_sim(alleles, mods, sex_det, num_gens, intro, efficiency, recomb_dist, \
                       fitness_cost, sterility_cost, add_intros, run_label, k, num_offspring, \
                       growth_factor, cross_dict, num_runs, file_name)

if __name__ == "__main__":
    main()

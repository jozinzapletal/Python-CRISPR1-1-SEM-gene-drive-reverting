# This module runs the CRISPR one-target location to one-target location using the equations generated and parameters
# set in the input file
#
# Model developed by Josef Zapletal (jozinzapletal@gmail.com)

import pandas as pd
import numpy as np
import math
import time
import CRISPR1_1_Equations


def CRISPR_1_1_Progression():

    cols = generate_allele_combinations(['w', 'v', 'u', 'r', 'g', 's'])
    results_temp = []

    # Reading and assigning values from Excel document
    M = pd.read_excel('CRISPR1-1 Input Parameters.xlsx')

    sigma = M.iloc[5, 3]
    lamda = M.iloc[6, 3]

    q2 = M.iloc[8, 3]
    q1 = M.iloc[9, 3]
    P2 = M.iloc[10, 3]
    P1 = M.iloc[11, 3]
    delta2 = M.iloc[12, 3]
    delta1 = M.iloc[13, 3]
    mortRateJuven = np.array([M.iloc[14, 3]] * 2 * len(cols))
    mortRateAdult1 = np.array([M.iloc[15, 3]] * 2 * len(cols))
    mortRateAdult2 = np.array([0.1] * 2 * len(cols))
    interval2 = 10
    mortRateAdult3 = np.array([0.1] * 2 * len(cols))
    interval3 = 20
    developTime = M.iloc[16, 3]
    gens = M.iloc[17, 3]
    span = gens

    male_alphas_init = M.iloc[0, 3:-1].values
    male_alphas = [i for i in male_alphas_init if not math.isnan(i)]  # refers to the male alpha values entered
    male_gammas_init = M.iloc[1, 3:-1].values
    male_gammas = [i for i in male_gammas_init if not math.isnan(i)]  # refers to the male gamma values entered

    female_alphas_init = M.iloc[2, 3:-1].values
    female_alphas = [i for i in female_alphas_init if not math.isnan(i)]  # refers to the male alpha values entered
    female_gammas_init = M.iloc[3, 3:-1].values
    female_gammas = [i for i in female_gammas_init if not math.isnan(i)]  # refers to the male gamma values entered

    # Assignment of initial conditions
    g0 = np.array(M.iloc[20:62, 3].values)

    # Assignment of initial fitness costs
    fitCost = np.array(M.iloc[20:62, 9].values)

    # recalculate fitness costs
    for index, val in enumerate(fitCost):
        mortRateAdult1[index] = mortRateAdult1[index] * (1 + fitCost[index])
        mortRateAdult2[index] = mortRateAdult2[index] * (1 + fitCost[index])
        mortRateAdult3[index] = mortRateAdult3[index] * (1 + fitCost[index])

        mortRateJuven[index] = mortRateJuven[index] * (1 + fitCost[index])

    # juvenile mortRate is calculated once, but encompasses the entire development time
    juvenile_mort_rate = [(1 - ((1 - i) ** developTime)) for i in mortRateJuven]
    juvenile_mort_rate = mortRateJuven


    # calculate the time period which to track carry-out adults
    min_adult_mort1 = min(mortRateAdult1)
    min_adult_mort2 = min(mortRateAdult2)
    min_adult_mort3 = min(mortRateAdult3)
    min_prop_surviving = 0.01
    prop_surviving = 1
    previous_days = 0
    while prop_surviving > min_prop_surviving:
        if previous_days < interval2:
            prop_surviving = prop_surviving * (1 - min_adult_mort1)
        elif interval2 <= previous_days <= interval3:
            prop_surviving = prop_surviving * (1 - min_adult_mort2)
        else:
            prop_surviving = prop_surviving * (1 - min_adult_mort3)
        previous_days += 1

    # RUN THE MODEL

    # indexes are day
    next_Generation = [[0 for x in range(2 * len(cols))] for y in range(span)]
    total_population_by_genotype = [[0 for x in range(2 * len(cols))] for y in range(span)]

    for alpha1 in male_alphas:
        for alpha2 in female_alphas:
            for gamma1 in male_gammas:
                for gamma2 in female_gammas:

                    beta1 = 1 - (alpha1 + gamma1)
                    beta2 = 1 - (alpha2 + gamma2)
                    epsilon1 = 0
                    epsilon2 = 0

                    # loop through each time step
                    start_time = time.time()
                    for T in range(span):
                        if T == 0:
                            proportionPop = convert_to_proportion(g0)
                            next_Generation[0] = CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost,
                                                 q1, q2, P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)
                        elif 0 < T:
                            proportionPop = convert_to_proportion(next_Generation[T-1])
                            next_Generation[T] = CRISPR1_1_Equations.CRISPR1_1_Rates(proportionPop, sigma, lamda, delta1, delta2, fitCost,
                                                 q1, q2, P1, P2, alpha1, alpha2, beta1, beta2, gamma1, gamma2, epsilon1, epsilon2)

                        # Store the values of the model into a list
                        total_population_by_genotype[T] = next_Generation[T]
                        total_population = sum(total_population_by_genotype[T])

                        male_w_count = 2*total_population_by_genotype[T][0] + sum(total_population_by_genotype[T][1:6])
                        male_v_count = total_population_by_genotype[T][1] + 2*total_population_by_genotype[T][6] + sum(total_population_by_genotype[T][7:11])
                        male_u_count = total_population_by_genotype[T][2] + total_population_by_genotype[T][7] + 2*total_population_by_genotype[T][11] + sum(total_population_by_genotype[T][12:15])
                        male_r_count = total_population_by_genotype[T][3] + total_population_by_genotype[T][8] + total_population_by_genotype[T][12] + 2*total_population_by_genotype[T][15] + sum(total_population_by_genotype[T][16:18])
                        male_g_count = total_population_by_genotype[T][4] + total_population_by_genotype[T][9] + total_population_by_genotype[T][13] + total_population_by_genotype[T][16] + 2*total_population_by_genotype[T][18] + total_population_by_genotype[T][19]
                        male_s_count = total_population_by_genotype[T][5] + total_population_by_genotype[T][10] + total_population_by_genotype[T][14] + total_population_by_genotype[T][17] + total_population_by_genotype[T][19] + 2*total_population_by_genotype[T][20]

                        female_w_count = 2*total_population_by_genotype[T][21] + sum(total_population_by_genotype[T][22:27])
                        female_v_count = total_population_by_genotype[T][22] + 2*total_population_by_genotype[T][27] + sum(total_population_by_genotype[T][28:32])
                        female_u_count = total_population_by_genotype[T][23] + total_population_by_genotype[T][28] + 2*total_population_by_genotype[T][32] + sum(total_population_by_genotype[T][33:36])
                        female_r_count = total_population_by_genotype[T][24] + total_population_by_genotype[T][29] + total_population_by_genotype[T][33] + 2*total_population_by_genotype[T][36] + sum(total_population_by_genotype[T][37:39])
                        female_g_count = total_population_by_genotype[T][25] + total_population_by_genotype[T][30] + total_population_by_genotype[T][34] + total_population_by_genotype[T][37] + 2*total_population_by_genotype[T][39] + total_population_by_genotype[T][40]
                        female_s_count = total_population_by_genotype[T][26] + total_population_by_genotype[T][31] + total_population_by_genotype[T][35] + total_population_by_genotype[T][38] + total_population_by_genotype[T][40] + 2*total_population_by_genotype[T][41]

                        w_count = (male_w_count + female_w_count) / (2*total_population)
                        v_count = (male_v_count + female_v_count) / (2*total_population)
                        u_count = (male_u_count + female_u_count) / (2*total_population)
                        r_count = (male_r_count + female_r_count) / (2*total_population)
                        g_count = (male_g_count + female_g_count) / (2*total_population)
                        s_count = (male_s_count + female_s_count) / (2*total_population)

                        output_line = [alpha1, gamma1, alpha2, gamma2, T, w_count, v_count, u_count, r_count, g_count, s_count, total_population]
                        results_temp.append(output_line)

                    print('Run complete in', time.time() - start_time, 'seconds.')

    # START HERE TO OUTPUT THE DATA
    cols = ['male_alpha', 'male_gamma', 'female_alpha', 'female_gamma', 'time', 'w', 'v', 'u', 'r', 'g', 's', 'total_population']
    results_df = pd.DataFrame(results_temp, columns=cols)

    results_df.to_excel('CRISPR1_1 Results.xlsx', index=False)


# converts the list of genotype counts to a proportion by
def convert_to_proportion(K):

    K = [float(i) for i in K]

    halfway_index = len(K)//2
    males = K[:halfway_index]
    females = K[halfway_index:]

    male_total = sum(males)

    if male_total != 0:
        males = [m / male_total for m in males]

    output_list = []
    output_list.extend(males)
    output_list.extend(females)

    return output_list


# generates a list of genotypes based on a list of alleles
def generate_allele_combinations(allele_list):

    # create all the combinations of alleles possible
    allele_1_sequence = []
    allele_2_sequence = []

    allele_dels = allele_list.copy()
    for allele in allele_list:
        for remaining_allele in allele_dels:
            allele_1_sequence.append(allele)
        allele_dels.remove(allele)

    allele_dels = allele_list.copy()
    for allele in allele_list:
        for remaining_allele in allele_dels:
            allele_2_sequence.append(remaining_allele)
        allele_dels.remove(allele)

    # combine the two alleles together into a single genotype
    two_allele_combinations = [i + j for i, j in zip(allele_1_sequence, allele_2_sequence)]

    return two_allele_combinations
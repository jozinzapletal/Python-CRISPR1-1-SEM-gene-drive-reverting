# This module generates the change equations that are required to calculate the changes from generation to generation

import os.path
from os import path

def CRISPR1_1_Eq_gen():

    #create all the possible two allele combinations
    allele_list = ['w', 'v', 'u', 'r', 'g', 's']
    two_allele_combinations = generate_allele_combinations(allele_list)

    # dictionary where all the equations will be stored
    equations = {}

    # loop through and assess the probabilities for the creation of alleles
    for result in two_allele_combinations:
        # count the types of alleles within the result genotype for easier comparison to the offspring genotype
        resultWcount = result.count('w')
        resultVcount = result.count('v')
        resultUcount = result.count('u')
        resultRcount = result.count('r')
        resultGcount = result.count('g')
        resultScount = result.count('s')

        equations[result] = ''

        for male in two_allele_combinations:
            # count the types of alleles within the male genotype
            maleWcount = male.count('w')
            maleVcount = male.count('v')
            maleUcount = male.count('u')
            maleRcount = male.count('r')
            maleGcount = male.count('g')
            maleScount = male.count('s')
            prob_gene_male = {}

            if (maleWcount == 1) and (maleGcount == 1):
                prob_gene_male['w'] = {'prob': '((alpha1/2)*(1-epsilon1)*(1+q1*P1)+(1-q1))'}
                prob_gene_male['v'] = {'prob': '((alpha1/2)*(epsilon1)*(1+q1*P1))'}
                prob_gene_male['u'] = {'prob': '((delta1/2)*q1*(1-P1))'}
                prob_gene_male['r'] = {'prob': '(((1-delta1)/2)*q1*(1-P1))'}
                prob_gene_male['g'] = {'prob': '((beta1/2)*(1+q1*P1))'}
                prob_gene_male['s'] = {'prob': '((gamma1/2)*(1+q1*P1))'}
            elif (maleWcount == 1) and (maleScount == 1):
                prob_gene_male['w'] = {'prob': '((1-q1)/2)'}
                prob_gene_male['u'] = {'prob': '((q1*(1-P1)*delta1)/2)'}
                prob_gene_male['r'] = {'prob': '((q1*(1-P1)*(1-delta1))/2)'}
                prob_gene_male['s'] = {'prob': '((1+q1*P1)/2)'}
            elif (maleVcount == 1) and (maleGcount == 1):
                prob_gene_male['w'] = {'prob': '((alpha1/2)*(1-epsilon1))'}
                prob_gene_male['v'] = {'prob': '((alpha1*epsilon1+1)/2)'}
                prob_gene_male['g'] = {'prob': '(beta1/2)'}
                prob_gene_male['s'] = {'prob': '(gamma1/2)'}
            elif (maleUcount == 1) and (maleGcount == 1):
                prob_gene_male['w'] = {'prob': '((alpha1/2)*(1-epsilon1))'}
                prob_gene_male['v'] = {'prob': '((alpha1/2)*(epsilon1))'}
                prob_gene_male['u'] = {'prob': '(1/2)'}
                prob_gene_male['g'] = {'prob': '(beta1/2)'}
                prob_gene_male['s'] = {'prob': '(gamma1/2)'}
            elif (maleRcount == 1) and (maleGcount == 1):
                prob_gene_male['w'] = {'prob': '((alpha1/2)*(1-epsilon1))'}
                prob_gene_male['v'] = {'prob': '((alpha1/2)*(epsilon1))'}
                prob_gene_male['r'] = {'prob': '(1/2)'}
                prob_gene_male['g'] = {'prob': '(beta1/2)'}
                prob_gene_male['s'] = {'prob': '(gamma1/2)'}
            elif (maleGcount == 1) and (maleScount == 1):
                prob_gene_male['w'] = {'prob': '((alpha1/2)*(1-epsilon1))'}
                prob_gene_male['v'] = {'prob': '((alpha1/2)*(epsilon1))'}
                prob_gene_male['g'] = {'prob': '(beta1/2)'}
                prob_gene_male['s'] = {'prob': '((gamma1+1)/2)'}
            elif maleGcount == 2:
                prob_gene_male['w'] = {'prob': '(alpha1*(1-epsilon1))'}
                prob_gene_male['v'] = {'prob': '(alpha1*epsilon1)'}
                prob_gene_male['g'] = {'prob': '(beta1)'}
                prob_gene_male['s'] = {'prob': '(gamma1)'}
            elif male[0] == male[1]:
                prob_gene_male[male[0]] = {'prob': '1'}
            else:
                prob_gene_male[male[0]] = {'prob': '1/2'}
                prob_gene_male[male[1]] = {'prob': '1/2'}

            for female in two_allele_combinations:
                # count the types of alleles within the female genotype
                femaleWcount = female.count('w')
                femaleVcount = female.count('v')
                femaleUcount = female.count('u')
                femaleRcount = female.count('r')
                femaleGcount = female.count('g')
                femaleScount = female.count('s')
                prob_gene_female = {}

                if (femaleWcount == 1) and (femaleGcount == 1):
                    prob_gene_female['w'] = {'prob': '((alpha2/2)*(1-epsilon2)*(1+q2*P2)+(1-q2))'}
                    prob_gene_female['v'] = {'prob': '((alpha2/2)*(epsilon2)*(1+q2*P2))'}
                    prob_gene_female['u'] = {'prob': '((delta2/2)*q2*(1-P2))'}
                    prob_gene_female['r'] = {'prob': '(((1-delta2)/2)*q2*(1-P2))'}
                    prob_gene_female['g'] = {'prob': '((beta2/2)*(1+q2*P2))'}
                    prob_gene_female['s'] = {'prob': '((gamma2/2)*(1+q2*P2))'}
                elif (femaleWcount == 1) and (femaleScount == 1):
                    prob_gene_female['w'] = {'prob': '((1-q2)/2)'}
                    prob_gene_female['u'] = {'prob': '((q2*(1-P2)*delta2)/2)'}
                    prob_gene_female['r'] = {'prob': '((q2*(1-P2)*(1-delta2))/2)'}
                    prob_gene_female['s'] = {'prob': '((1+q2*P2)/2)'}
                elif (femaleVcount == 1) and (femaleGcount == 1):
                    prob_gene_female['w'] = {'prob': '((alpha2/2)*(1-epsilon2))'}
                    prob_gene_female['v'] = {'prob': '((alpha2*epsilon2+1)/2)'}
                    prob_gene_female['g'] = {'prob': '(beta2/2)'}
                    prob_gene_female['s'] = {'prob': '(gamma2/2)'}
                elif (femaleUcount == 1) and (femaleGcount == 1):
                    prob_gene_female['w'] = {'prob': '((alpha2/2)*(1-epsilon2))'}
                    prob_gene_female['v'] = {'prob': '((alpha2/2)*(epsilon2))'}
                    prob_gene_female['u'] = {'prob': '(1/2)'}
                    prob_gene_female['g'] = {'prob': '(beta2/2)'}
                    prob_gene_female['s'] = {'prob': '(gamma2/2)'}
                elif (femaleRcount == 1) and (femaleGcount == 1):
                    prob_gene_female['w'] = {'prob': '((alpha2/2)*(1-epsilon2))'}
                    prob_gene_female['v'] = {'prob': '((alpha2/2)*(epsilon2))'}
                    prob_gene_female['r'] = {'prob': '(1/2)'}
                    prob_gene_female['g'] = {'prob': '(beta2/2)'}
                    prob_gene_female['s'] = {'prob': '(gamma2/2)'}
                elif (femaleGcount == 1) and (femaleScount == 1):
                    prob_gene_female['w'] = {'prob': '((alpha2/2)*(1-epsilon2))'}
                    prob_gene_female['v'] = {'prob': '((alpha2/2)*(epsilon2))'}
                    prob_gene_female['g'] = {'prob': '(beta2/2)'}
                    prob_gene_female['s'] = {'prob': '((gamma2+1)/2)'}
                elif femaleGcount == 2:
                    prob_gene_female['w'] = {'prob': '(alpha2*(1-epsilon2))'}
                    prob_gene_female['v'] = {'prob': '(alpha2*epsilon2)'}
                    prob_gene_female['g'] = {'prob': '(beta2)'}
                    prob_gene_female['s'] = {'prob': '(gamma2)'}
                elif female[0] == female[1]:
                    prob_gene_female[female[0]] = {'prob': '1'}
                else:
                    prob_gene_female[female[0]] = {'prob': '1/2'}
                    prob_gene_female[female[1]] = {'prob': '1/2'}

                # COMBINE INTO OFFSPRING GENOTYPE

                prob_gene_offspring = {}
                for male_allele in prob_gene_male.keys():
                    for female_allele in prob_gene_female.keys():
                        offspring = male_allele + female_allele
                        # count the types of alleles within the offspring genotype
                        offspringWcount = offspring.count('w')
                        offspringVcount = offspring.count('v')
                        offspringUcount = offspring.count('u')
                        offspringRcount = offspring.count('r')
                        offspringGcount = offspring.count('g')
                        offspringScount = offspring.count('s')

                        if offspringWcount == resultWcount and offspringVcount == resultVcount \
                                and offspringUcount == resultUcount and offspringRcount == resultRcount \
                                and offspringGcount == resultGcount and offspringScount == resultScount:

                            prob_gene_offspring = prob_gene_male[male_allele]['prob'] + '*' + \
                                                    prob_gene_female[female_allele]['prob']

                            combined_prob_funct = 'J[' + str(two_allele_combinations.index(male)) + ']*' + prob_gene_offspring + \
                                                  '*J[' + str(two_allele_combinations.index(female) + len(two_allele_combinations)) + ']'

                            # store the resulting equation into the final dictionary
                            if equations[result] == '':
                                equations[result] = combined_prob_funct
                            else:
                                equations[result] = equations[result] + ' + ' + combined_prob_funct

        # add in the fitness cost, sex ratios, and laying rate

    # check if the file exists - create new file if false, replace if true
    file_name = 'CRISPR1_1_Equations.py'
    if path.exists(file_name):
        os.remove(file_name)
        f = open(file_name, 'w')
        f.write('def CRISPR1_1_Rates(J, sigma, lamda, delta1, delta2, fitCost, q1, q2, P1, P2, alpha1, alpha2, beta1, '
                'beta2, gamma1, gamma2, epsilon1, epsilon2):\n'
                '    rates = []\n\n')
        eq_num = 0
        for eq in equations:
            f.write('    rates.append((1-fitCost[' + str(eq_num) + '])*sigma*lamda*(' + equations[eq] + '))\n')
            eq_num += 1
        for eq in equations:
            f.write('    rates.append((1-fitCost[' + str(eq_num) + '])*(1-sigma)*lamda*(' + equations[eq] + '))\n')
            eq_num += 1
        f.write('\n    return rates\n')
        f.close()

    else:
        f = open(file_name, 'w')
        f.write('def CRISPR1_1_Rates(J, sigma, lamda, delta1, delta2, fitCost, q1, q2, P1, P2, alpha1, alpha2, beta1, '
                'beta2, gamma1, gamma2, epsilon1, epsilon2):\n'
                '    rates = []\n\n')
        eq_num = 0
        for eq in equations:
            f.write('    rates.append((1-fitCost[' + str(eq_num) + '])*sigma*lamda*(' + equations[
                eq] + '))\n')
            eq_num += 1
        for eq in equations:
            f.write('    rates.append((1-fitCost[' + str(eq_num) + '])*(1-sigma)*lamda*(' + equations[
                eq] + '))\n')
            eq_num += 1
        f.write('\n    return rates\n\n')
        f.close()

        print('Equation generation complete.')


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

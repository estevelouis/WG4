from operator import *
import math
from math import inf
import scipy
import numpy as np
import copy
import pandas as pd
import array
import wasserstein
import sys
import matplotlib.pyplot as plt
import copy
from python_tsp.exact import solve_tsp_dynamic_programming
from ups_lingdiv_utils.zipfian_series import zipfian_series
if "--no-gensim" not in sys.argv:
    import gensim.downloader
    glove_vectors = gensim.downloader.load('glove-wiki-gigaword-50')

LOG_BASE = 2.0

ALLOW_RATIO = True
ALLOW_RATIO_IF_MATH_INF = False

CONVERT_PROPORTION_TO_ABUNDANCE = True
CONVERT_PROPORTION_TO_ABUNDANCE_RATIO = 1000

ENABLE_MAX_DIMENSIONALITY = True
MAX_DIMENSIONALITY = 2

ROUND_TO_CORRECT_FP_INACCURACY_AXIOM1 = True
ROUND_TO_CORRECT_FP_INACCURACY_AXIOM2 = True
ROUND_TO_CORRECT_FP_INACCURACY_AXIOM3 = True
ROUND_TO_CORRECT_FP_INACCURACY_AXIOM4 = True
ROUND_TO_CORRECT_FP_INACCURACY_AXIOM5 = True
ROUND_TO_CORRECT_FP_INACCURACY_D1 = True
ROUND_TO_CORRECT_FP_INACCURACY_D2 = True
ROUND_TO_CORRECT_FP_INACCURACY_D3 = True
ROUND_TO_CORRECT_FP_INACCURACY_D4 = True
ROUND_TO_CORRECT_FP_INACCURACY_P1 = True
ROUND_TO_CORRECT_FP_INACCURACY_P2 = True
ROUND_TO_CORRECT_FP_INACCURACY_SW2 = True
ROUND_TO_CORRECT_FP_INACCURACY_SW3 = True
ROUND_TO_CORRECT_FP_INACCURACY_SW6 = True
ROUND_TO_CORRECT_FP_INACCURACY_SW7 = True
ROUND_TO_CORRECT_FP_INACCURACY_SW13 = True
ROUND_TO_CORRECT_FP_INACCURACY_SW14 = True
ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS = 8

def to_ratio(a):
    if (not ALLOW_RATIO):
        return array.array("d", a)
    if (not ALLOW_RATIO_IF_MATH_INF) and (math.inf in a or (-math.inf) in a):
        return array.array("d", a)
    sum_ = sum(a)
    return array.array("d", [i/sum_ for i in a])

def to_abundance(a):
    return [i * CONVERT_PROPORTION_TO_ABUNDANCE_RATIO for i in a]

def get_vector_with_ratio(key_, ratio=1.0):
    vec = glove_vectors[key_]
    if ENABLE_MAX_DIMENSIONALITY:
        vec = vec[:MAX_DIMENSIONALITY]
    vec = vec * ratio
    return vec

def df_species_count(a):
    a = to_ratio(a)
    result = (len(a) - 1)
    assert(type(result) == int or type(result) == float)
    return result

def df_shannon_entropy(a):
    a = to_ratio(a)
    result = 0.0
    for i in a:
        if i == 0.0:
            continue
        result += i * math.log(i, LOG_BASE)
    result = - result
    assert(type(result) == int or type(result) == float)
    return result

def df_simpson_index(a):
    a = to_ratio(a)
    result = 1.0
    sum_ = 0.0
    for i in a:
        sum_ += i ** 2
    result -= sum_
    assert(type(result) == int or type(result) == float)
    return result

def df_richness(a):
    a = to_ratio(a)
    result = len(a)
    return result

def df_richness_ht0(a):
    a = to_ratio(a)
    result = 0
    for i in a:
        if i > 0:
            result += 1
    return result

def hill_number(a, base=2):
    a = to_ratio(a)
    if base == 1:
        result = df_shannon_entropy(a)
    else:
        result = 0.0
        for i in a:
            result += i ** base
        result = result ** (1 / (1-base))
    return result


def df_hill_number_base_two(a):
    a = to_ratio(a)
    return hill_number(a, base=2)

def df_hill_number_base_e(a):
    a = to_ratio(a)
    return hill_number(a, base=math.e)

def df_hill_number_base_ten(a):
    a = to_ratio(a)
    return hill_number(a, base=10)

def hill_evenness_x_y(a, x=1, y=0):
    a = to_ratio(a)
    upper_half = hill_number(a, base=x)
    lower_half = hill_number(a, base=y)
    result = upper_half / lower_half
    return result

def df_hill_evenness_1_0(a):
    a = to_ratio(a)
    return hill_evenness_x_y(a, x=1, y=0)

def df_hill_evenness_2_1(a):
    a = to_ratio(a)
    return hill_evenness_x_y(a, x=2, y=1)

def df_morales_et_al_D1_shannon_diversity(a):
    a = to_ratio(a)
    result = 1.0
    for i in a:
        if i == 0.0:
            continue
        result *= (i ** i)
    result = 1.0 / result
    return result
    
def df_morales_et_al_D2_herfindahl_diversity(a):
    a = to_ratio(a)
    result = 0.0
    for i in a:
        if i == 0.0:
            continue
        result += (i ** 2)
    result = 1.0 / result
    return result

def df_morales_et_al_Dinf_berger_diversity(a):
    a = to_ratio(a)
    max_ = max(a)
    result = 1.0 / max_
    return result
    
def df_herfindahl_hirschman_index(a):
    a = to_ratio(a)
    result = 0.0
    for i in a:
        if i == 0.0:
            continue
        result += (i ** 2)
    return result

def df_berger_parker_index(a):
    a = to_ratio(a)
    result = max(a)
    return result

def df_ineq_gini_coefficient_old(a):
    a = to_ratio(a)
    result = 0.0
    for i in range(len(a)):
        for j in range(len(a)):
            result += abs(a[i] - a[j])
    result *= (1.0 / (2 * len(a)))
    return result

def df_hr_ineq_gini_coefficient_new(a):
    # see Hurley and Rickard (2009)
    a = to_ratio(a)
    a = sorted(a) # this measure requires positively increasing values
    norm_1 = sum(a)
    sum_ = sum((a[k] / norm_1) * ((len(a) - k + (1.0/2.0)) / len(a)) for k in range(len(a)))
    result = 1.0 - 2.0 * sum_
    return result

def multilingual_df_hr_ineq_gini_coefficient_new(a, consideration_for_each=[], norm_1=1.0, norm_2=1.0):
    print(len(a), len(consideration_for_each))
    if len(consideration_for_each) != len(a):
        print("warning: different length between a and consideration for each")
        consideration_for_each = [1.0] * len(a)
    #a = to_ratio(a) # doing it custom instead
    sum_ = sum(a)
    a = [i / sum_ for i in a]
    a = [(i, j) for i, j in zip(a, consideration_for_each)]
    a = sorted(a, key=lambda x:x[0]) # this measure requires positively increasing values
    for i in a:
        print(i)
    sum_ = sum((a[k][0] / norm_1) * (a[k][1] / norm_2) * ((len(a) - k + (1.0/2.0)) / len(a)) for k in range(len(a)))
    result = 1.0 - 2.0 * sum_
    return result

def multilingual_df_hr_ineq_gini_coefficient_new_v2(a, consideration_for_each=[], k_limit=1, *args):
    print(len(a), len(consideration_for_each))
    print("k_limit", k_limit)
    if len(consideration_for_each) != len(a):
        print("warning: different length between a and consideration for each")
        consideration_for_each = [1.0] * len(a)
    a = [i * j for i, j in zip(a, consideration_for_each)]
    norm_1 = sum(a)
    a = sorted(a) # this measure requires positively increasing values
    for i in a:
        print(i)
    sum_ = sum((a[k] / norm_1) * ((len(a) - k + (1.0/2.0)) / len(a)) for k in range(k_limit))
    result = 1.0 - 2.0 * sum_
    return result

def df_hr_l0(a):
    return len(list(filter(lambda x:x==0.0, a)))

def _df_hr_l_epsilon(a, epsilon):
    return len(list(filter(lambda x:x<=epsilon, a)))

def df_hr_l_epsilon1(a):
    return _df_hr_l_epsilon(a, epsilon=1.0)
def df_hr_l_epsilon2(a):
    return _df_hr_l_epsilon(a, epsilon=2.0)
def df_hr_l_epsilon3(a):
    return _df_hr_l_epsilon(a, epsilon=3.0)
def df_hr_l_epsilon4(a):
    return _df_hr_l_epsilon(a, epsilon=4.0)
def df_hr_l_epsilon5(a):
    return _df_hr_l_epsilon(a, epsilon=5.0)

def df_hr_negative_l1(a):
    a = to_ratio(a)
    sum_ = sum(a)
    return -sum_

def _df_hr_negative_lp(a, p):
    assert(0.0 < p < 1.0)
    a = to_ratio(a)
    sum_ = sum(i ** p for i in a)
    return -(sum_ ** (1.0 / p))

def df_hr_negative_l010(a):
    return _df_hr_negative_lp(a, 0.10)
def df_hr_negative_l050(a):
    return _df_hr_negative_lp(a, 0.50)
def df_hr_negative_l090(a):
    return _df_hr_negative_lp(a, 0.90)

def df_hr_l2_over_l1(a):
    a = to_ratio(a)
    upper_sum = sum(i ** 2.0 for i in a)
    lower_sum = sum(a)
    return math.sqrt(upper_sum) / lower_sum

def _df_hr_negative_tanh_a_b(a, alpha, beta):
    a = to_ratio(a)
    sum_ = sum(math.tanh((alpha*i) ** beta) for i in a)
    return -sum_

def df_hr_negative_tanh_a2_b2(a):
    return _df_hr_negative_tanh_a_b(a, alpha=2.0, beta=2.0)
def df_hr_negative_tanh_a2_b3(a):
    return _df_hr_negative_tanh_a_b(a, alpha=2.0, beta=3.0)
def df_hr_negative_tanh_a2_b4(a):
    return _df_hr_negative_tanh_a_b(a, alpha=2.0, beta=4.0)
def df_hr_negative_tanh_a3_b2(a):
    return _df_hr_negative_tanh_a_b(a, alpha=3.0, beta=2.0)
def df_hr_negative_tanh_a3_b3(a):
    return _df_hr_negative_tanh_a_b(a, alpha=3.0, beta=3.0)
def df_hr_negative_tanh_a3_b4(a):
    return _df_hr_negative_tanh_a_b(a, alpha=3.0, beta=4.0)
def df_hr_negative_tanh_a4_b2(a):
    return _df_hr_negative_tanh_a_b(a, alpha=4.0, beta=2.0)
def df_hr_negative_tanh_a4_b3(a):
    return _df_hr_negative_tanh_a_b(a, alpha=4.0, beta=3.0)
def df_hr_negative_tanh_a4_b4(a):
    return _df_hr_negative_tanh_a_b(a, alpha=4.0, beta=4.0)

def df_hr_negative_log(a):
    a = to_ratio(a)
    sum_ = sum(math.log(1.0 + i, LOG_BASE) for i in a)
    return -sum_

def df_hr_k4(a):
    a = to_ratio(a)
    upper_sum = sum(i ** 4.0 for i in a)
    lower_sum = sum(i ** 2.0 for i in a)
    result = upper_sum / (lower_sum ** 2.0)
    return result

def _df_hr_negative_lp_no_negative(a, p):
    assert(p < 0.0)
    a = to_ratio(a)
    sum_ = sum(i ** p for i in filter(lambda x:x!=0, a))
    return -sum_

def df_hr_negative_lneg010_no_negative(a):
    return _df_hr_negative_lp_no_negative(a, p=-0.10)
def df_hr_negative_lneg050_no_negative(a):
    return _df_hr_negative_lp_no_negative(a, p=-0.50)
def df_hr_negative_lneg090_no_negative(a):
    return _df_hr_negative_lp_no_negative(a, p=-0.90)

def df_hr_hg(a):
    a = to_ratio(a)
    sum_ = sum(math.log(i ** 2.0, LOG_BASE) for i in a if i != 0.0)
    return -sum_

def df_hr_h_prime_s(a):
    a = to_ratio(a)
    sum_ = sum(i * math.log(i ** 2.0, LOG_BASE) for i in a if i != 0.0)
    return -sum_

def df_hr_hoyer(a):
    a = to_ratio(a)
    upper_sum = sum(a)
    lower_sum = sum(i ** 2.0 for i in a)
    sqrt_n = math.sqrt(len(a))
    result = (sqrt_n - (upper_sum / math.sqrt(lower_sum))) * ((sqrt_n - 1.0) ** -1.0)
    return result

def df_shannon_evenness(a):
    # uncertain, but based on description from Morales et al. (2020), should be like this?
    a = to_ratio(a)

    observed_entropy = df_shannon_entropy(a)
    a_ = [1.0/len(a) for i in a]
    maximum_entropy = df_shannon_entropy(a_)
    result = observed_entropy / maximum_entropy
    return result

def concentration_ratio(a, n=1):
    # original equation in base-1 does [1,n], but as here we are in base-0 we do [0,n[
    a = to_ratio(a)
    b = sorted(a)
    result = 0.0
    for i in range(min(n, len(b))):
        result += b[i]
    return result

def df_concentration_ratio_n1(a):
    a = to_ratio(a)
    return concentration_ratio(a, n=1)
    
def df_concentration_ratio_n2(a):
    a = to_ratio(a)
    return concentration_ratio(a, n=2)
    
"""
def df_junge1994_page20(a):
    a = to_ratio(a)
    try:
        n = len(a)
        sum_ = sum((i ** 2) - 1 for i in a)
        right = math.sqrt(n - 1) - math.sqrt(n * sum_)
        left = (1 / math.sqrt(n))
        return left * right
    except Exception as e:
        print(e)
        return None
"""

def df_junge1994_page22(a):
    a = to_ratio(a) # fix
    sum_ = 0.0
    for i in range(len(a)):
        sum_ += a[i] ** 2
    result = 1 - (sum_ ** (1/2))
    return result

def df_brillouin_div(a):
    # no normalisation
    # see Laxton (1978)
    n = len(a)

    """
    try:
        upper = math.factorial(n)
        lower = 1
        for i in a:
            lower *= math.factorial(int(i))
        print(upper, lower)
        result = math.log(upper / lower, LOG_BASE)
        return result
    except Exception as e:
        print(e)
        return None
    """

    sum_left = 0.0
    sum_right = 0.0
    for i in range(n):
        sum_left += math.log(i+1, LOG_BASE)
        sum_right += math.log(math.factorial(int(a[i])), LOG_BASE)
    return sum_left - sum_right


def df_mcintosh_index(a):
    # see Laxton (1978)
    a = to_ratio(a)
    sum_ = sum(i ** 2 for i in a)
    result = 1 - math.sqrt(sum_)
    return result

def df_sw_simpson_dominance_index(a):
    a = to_ratio(a)
    sum_ = 0.0
    for i in a:
        sum_ += i ** 2
    return sum_

# should len(a) only be non-zero species?
def df_sw_entropy_over_log_n_species_pielou1975(a):
    a = to_ratio(a)
    entropy_ = df_shannon_entropy(a)
    return entropy_ / math.log(len(a), math.e)

def renyi_entropy(a, alpha=1.0):
    a = to_ratio(a)
    if alpha == 1.0:
        sum_ = sum(i * math.log(i, LOG_BASE) for i in a)
        return -sum_
    else:
        sum_ = sum(i ** alpha for i in a)
        return (1.0 / (1.0 - alpha)) * math.log(sum_, LOG_BASE)
    
def df_sw_e_heip(a):
    a = to_ratio(a)
    entropy_ = df_shannon_entropy(a)
    result = ((math.e ** entropy_) - 1) / (len(a) - 1)
    return result

def _df_sw_e_heip_custom_scale(a, alpha):
    a = to_ratio(a)
    entropy_ = renyi_entropy(a, alpha)
    result = ((math.e ** entropy_) - 1) / (len(a) - 1)
    return result

def df_sw_e_heip_custom_scale_alpha_095(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.95)
def df_sw_e_heip_custom_scale_alpha_094(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.94)
def df_sw_e_heip_custom_scale_alpha_093(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.93)
def df_sw_e_heip_custom_scale_alpha_092(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.92)
def df_sw_e_heip_custom_scale_alpha_091(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.91)
def df_sw_e_heip_custom_scale_alpha_090(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.90)
def df_sw_e_heip_custom_scale_alpha_089(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.89)
def df_sw_e_heip_custom_scale_alpha_088(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.88)
def df_sw_e_heip_custom_scale_alpha_087(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.87)
def df_sw_e_heip_custom_scale_alpha_086(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.86)
def df_sw_e_heip_custom_scale_alpha_085(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.85)
def df_sw_e_heip_custom_scale_alpha_084(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.84)
def df_sw_e_heip_custom_scale_alpha_083(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.83)
def df_sw_e_heip_custom_scale_alpha_082(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.82)
def df_sw_e_heip_custom_scale_alpha_081(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.81)
def df_sw_e_heip_custom_scale_alpha_080(a):
    return _df_sw_e_heip_custom_scale(a, alpha=0.80)
    
def df_sw_e_one_minus_D(a):
    a = to_ratio(a)
    dominance = df_sw_simpson_dominance_index(a)
    result = (1 - dominance) / (1 - (1/len(a)))
    return result

def df_sw_e_one_over_D_williams1964(a):
    a = to_ratio(a)
    dominance = df_sw_simpson_dominance_index(a)
    result = (1/dominance) / len(a)
    return result

def df_sw_e_minus_ln_D_pielou1977(a):
    a = to_ratio(a)
    dominance = df_sw_simpson_dominance_index(a)
    result = (-math.log(dominance, math.e)) / math.log(len(a), math.e)
    return result

def df_sw_f_2_1_alatalo1981(a):
    a = to_ratio(a)
    dominance = df_sw_simpson_dominance_index(a)
    entropy_ = df_shannon_entropy(a)
    result = ((1/dominance) - 1) / ((math.e ** entropy_) - 1)
    return result

def df_sw_g_2_1_molinari1989(a):
    a = to_ratio(a)
    f_2_1 = df_sw_f_2_1_alatalo1981(a)
    if f_2_1 > math.sqrt(1/2):
        result = f_2_1 * 0.636611 * math.asin(f_2_1)
    else:
        result = f_2_1 ** 3
    return result

def df_sw_o_bulla1994(a):
    a = to_ratio(a)
    sum_ = 0.0
    for i in a:
        sum_ += min(i, 1/len(a))
    return sum_

def df_sw_e_bulla1994(a):
    a = to_ratio(a)
    O_ = df_sw_o_bulla1994(a)
    result = (O_ - (1/len(a))) / (1 - (1/len(a)))
    return result

def df_sw_e_mci_pielou1969(a):
    #a = to_ratio(a) # explicitly disabled here as we need the number of items
    sum_x = 0
    sum_x_square = 0
    for i in a:
        sum_x += i
        sum_x_square += i ** 2
    result = (sum_x - math.sqrt(sum_x_square)) / (sum_x - (sum_x / math.sqrt(len(a))))
    return result

def df_sw_e_prime_camargo1993(a):
    a = to_ratio(a)
    sum_ = 0.0
    n = len(a)
    for i in range(n):
        for j in range(i+1,n):
            sum_ = add(sum_, abs(a[i] - a[j]) / n)
    result = 1 - sum_
    return result

def df_sw_e_var_smith_and_wilson1996_original(a):
    #a = to_ratio(a) # explicitly disabled
    sum_ = 0.0
    n = len(a)

    inner_sum = sum(map(lambda x:math.log(x, math.e) / n if x != 0.0 and n != 0 else 0.0, a))
    outer_sum = sum(map(lambda x:((math.log(x, math.e) - inner_sum) ** 2) / n if x != 0.0 and n != 0 else 0.0, a))
    result = 1.0 - ((2.0/math.pi) * math.atan(outer_sum))

    return result

def distance_minkowski_p(a, b, p):
    assert(len(a) == len(b))
    sum_ = 0.0
    for i in range(len(a)):
        sum_ += abs(a[i] - b[i]) ** p
    return (sum_ ** (1/p))

def distance_minkowski_p1(a, b):
    return distance_minkowski_p(a, b, p=1.0)

def distance_minkowski_p2(a, b):
    return distance_minkowski_p(a, b, p=2.0)

def _intraset_disparity_ricotta_szeidl(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0, alpha=2.0):
    a = to_ratio(a)
    all_vectors = list(get_vector_with_ratio(i, vector_ratio) for i in keys_)
    all_vectors = np.array(all_vectors)

    upper_sum = 0.0
    for i in range(len(a)):
        inner_sum = 0.0
        for j in range(len(a)):
            if i == j:
                continue
            inner_sum += distance_function(all_vectors[i], all_vectors[j]) * a[j]
        upper_sum += a[i] * ((1 - inner_sum) ** (alpha - 1))

    q_alpha = (1 - upper_sum) / (alpha - 1)
    return q_alpha
    

def _intraset_disparity_weitzman(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0):
    def weitzman(matrix):
        path = []
        if __debug__:
            print(matrix)
        argmin_ = np.argmin(matrix) # yields argmin if flattened array 
        argmin_dim1_index = math.floor(argmin_ / len(keys_))
        argmin_dim2_index = argmin_ % len(keys_)

        min_value = matrix[argmin_dim1_index][argmin_dim2_index]
        if min_value == np.inf:
            return (0.0, path) # empty path
        else:
            option_a = copy.deepcopy(matrix)
            option_b = copy.deepcopy(matrix)
            
            option_a[:, argmin_dim1_index] = np.inf
            option_a[argmin_dim1_index, :] = np.inf
            option_b[:, argmin_dim2_index] = np.inf
            option_b[argmin_dim2_index, :] = np.inf

            option_a_res = weitzman(option_a)
            option_b_res = weitzman(option_b)
            if option_a_res[0] >= option_b_res[0]:
                path.append(argmin_dim1_index)
                chosen_result = option_a_res
            else:
                path.append(argmin_dim2_index)
                chosen_result = option_b_res
            path += chosen_result[1]
            return (min_value + chosen_result[0], path)

    all_vectors = list(get_vector_with_ratio(i, vector_ratio) for i in keys_)
    all_vectors = np.array(all_vectors)
    
    distance_matrix = np.ones((len(keys_), len(keys_))) * np.inf
    for i in range(len(keys_)):
        for j in range(i+1, len(keys_)):
            distance_matrix[i][j] = distance_function(all_vectors[i], all_vectors[j])

    return weitzman(distance_matrix)

"""
def _intraset_disparity_lexicographic(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0):
    # inaccurate, see rework in C for better implementation
    def lexicographic(matrix):
        path = []
        if __debug__:
            print(matrix)
        argmin_ = np.argmin(matrix) # yields argmin if flattened array 
        argmin_dim1_index = math.floor(argmin_ / len(keys_))
        argmin_dim2_index = argmin_ % len(keys_)

        min_value = matrix[argmin_dim1_index][argmin_dim2_index]
        if min_value == np.inf:
            return (0.0, path) # empty path
        else:
            option_a = copy.deepcopy(matrix)
            #option_b = copy.deepcopy(matrix)
            
            option_a[:, argmin_dim1_index] = np.inf
            option_a[argmin_dim1_index, :] = np.inf
            #option_b[:, argmin_dim2_index] = np.inf
            #option_b[argmin_dim2_index, :] = np.inf

            option_a_res = lexicographic(option_a)
            #option_b_res = lexicographic(option_b)
            chosen_result = option_a_res
            path += chosen_result[1]
            return (min_value + chosen_result[0], path)

    all_vectors = list(get_vector_with_ratio(i, vector_ratio) for i in keys_)
    all_vectors = np.array(all_vectors)
    
    distance_matrix = np.ones((len(keys_), len(keys_))) * np.inf
    for i in range(len(keys_)):
        for j in range(i+1, len(keys_)):
            distance_matrix[i][j] = distance_function(all_vectors[i], all_vectors[j])

    return lexicographic(distance_matrix)
"""

def _intraset_disparity_tsp(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0):
    all_vectors = list(get_vector_with_ratio(i, vector_ratio) for i in keys_)
    all_vectors = np.array(all_vectors)
    
    # modified, different from Weitzman
    distance_matrix = np.zeros((len(keys_), len(keys_)))
    for i in range(len(keys_)):
        for j in range(len(keys_)):
            distance_matrix[i][j] = distance_function(all_vectors[i], all_vectors[j])
            distance_matrix[j][i] = distance_function(all_vectors[j], all_vectors[i])

    permutation, distance = solve_tsp_dynamic_programming(distance_matrix)
    extra_distance = distance_function(all_vectors[permutation[0]], all_vectors[permutation[-1]])
    return (distance, permutation, extra_distance)


def _intraset_disparity_stirling2007_delta(a, keys_=[], alpha=1.0, beta=1.0, distance_function=distance_minkowski_p2, vector_ratio=1.0):
    assert(alpha == 0.0 or len(a) == len(keys_))
    a = to_ratio(a)

    a_powers = []

    def compute_powers():
        nonlocal a_powers
        a_powers = array.array("d", (i ** beta for i in a))

    if alpha == 0.0:
        if beta == 0.0:
            n = len(a)
            return (n * (n-1)) / 2
        else:
            compute_powers()
            sum_ = 0.0
            for i in range(len(a)):
                for j in range(i+1,len(a)):
                    sum_ = add(sum_, mul(a_powers[i], a_powers[j]))
            return sum_
    else:
        if beta == 0.0:
            #sum_ = scipy.spatial.distance.pdist(keys_, metric=lambda a_, b_:distance_function(get_vector_with_ratio(a_, vector_ratio), get_vector_with_ratio(b_, vector_ratio)) ** alpha)
            sum_ = 0.0
            for i in range(len(a)):
                for j in range(i+1, len(a)):
                    sum_ = add(sum_, distance_function(get_vector_with_ratio(keys_[i], vector_ratio), get_vector_with_ratio(keys_[j], vector_ratio)) ** alpha)
            return sum_
        else:
            sum_ = 0.0
            for i in range(len(a)):
                keys_i_found = keys_[i] in glove_vectors
                for j in range(i+1, len(a)):
                    distance_result = 1.0
                    if alpha != 0 and keys_i_found and keys_[j] in glove_vectors:
                        distance_result = distance_function(get_vector_with_ratio(keys_[i], vector_ratio), get_vector_with_ratio(keys_[j], vector_ratio)) ** alpha
                    balance_result = 1.0
                    if beta != 0:
                        balance_result = (a[i] * a[j]) ** beta
                    sum_ += distance_result * balance_result
            return sum_

def df_intraset_disparity_stirling2007_delta_alpha0_beta0(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0):
    return _intraset_disparity_stirling2007_delta(a, keys_, alpha=0.0, beta=0.0, distance_function=distance_function, vector_ratio=vector_ratio)

def df_intraset_disparity_stirling2007_delta_alpha0_beta1(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0):
    return _intraset_disparity_stirling2007_delta(a, keys_, alpha=0.0, beta=1.0, distance_function=distance_function, vector_ratio=vector_ratio)

def _df_intraset_disparity_stirling2007_delta_alpha1_beta0(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0):
    return _intraset_disparity_stirling2007_delta(a, keys_, alpha=1.0, beta=0.0, distance_function=distance_function, vector_ratio=vector_ratio)

def _df_intraset_disparity_stirling2007_delta_alpha1_beta1(a, keys_=[], distance_function=distance_minkowski_p2, vector_ratio=1.0):
    return _intraset_disparity_stirling2007_delta(a, keys_, alpha=1.0, beta=1.0, distance_function=distance_function, vector_ratio=vector_ratio)

def _intraset_disparity_leinster_and_cobbold2012(a, keys_=[], alpha=1.0, beta=1.0, distance_function=distance_minkowski_p2, vector_ratio=1.0):
    assert(alpha == 0.0 or len(a) == len(keys_))
    a = to_ratio(a)

    if alpha == 1.0:
        prod_ = 1.0
        for i in range(len(a)):
            Zp = 0.0
            for j in range(len(a)):
                Zp += distance_function(get_vector_with_ratio(keys_[i], vector_ratio), get_vector_with_ratio(keys_[j], vector_ratio)) * a[j]
            prod_ *= Zp ** a[j]
        result = 1.0 / prod_
        return result
    else:
        sum_ = 0.0
        for i in range(len(a)):
            Zp = 0.0
            for j in range(len(a)):
                Zp += distance_function(get_vector_with_ratio(keys_[i], vector_ratio), get_vector_with_ratio(keys_[j], vector_ratio)) * a[j]
            sum_ += Zp ** (alpha-1)
        result = sum_ ** (1/(1-alpha))
        return result

"""
def _intraset_disparity_chao_and_jost2014_unified_functional(a, keys_=[], alpha=1.0, beta=1.0, distance_function=distance_minkowski_p2, vector_ratio=1.0):
    a = to_ratio(a)
    dist_matrix = np.zeros((len(a), len(a)), dtype=np.float64)
    for i in range(len(a)):
        for j in range(len(a)):
            dist_matrix[i][j] = distance_function(get_vector_with_ratio(keys_[i], vector_ratio), get_vector_with_ratio(keys_[j], vector_ratio)) * a[i] * a[j]

    if alpha == 1.0:
        sum_ = 0.0
        for i in range(len(a)):
            for j in range(len(a)):
                ratio = (a[i] * a[j]) / sum_matrix
                sum_ += dist_matrix[i][j] * ratio * math.log(ratio)
        result = math.exp(-sum_)
        return result
    else:
        sum_ = 0.0
        sum_matrix = np.sum(dist_matrix)
        for i in range(len(a)):
            for j in range(len(a)):
                sum_ += dist_matrix[i][j] * (((a[i] * a[j]) / sum_matrix) ** alpha)
    
        result = sum_ ** (1/(1-alpha))
        return result
"""

def _intraset_disparity_chao_and_jost2014_unified_species(a, keys_=[], alpha=1.0, beta=1.0, distance_function=distance_minkowski_p2, vector_ratio=1.0):
    a = to_ratio(a)

    if alpha == 1.0:
        print("not implemented")
        pass
    else:
        sum_ = 0.0
        for i in range(len(a)):
            sum_ += a[i] ** alpha
    
        result = sum_ ** (1/(1-alpha))
        return result

def main():
    axiom1_samples = [
        {"a": [0.25, 0.75], "b": [0.75, 0.25]},
        {"a": [0.40, 0.60], "b": [0.60, 0.40]},
        {"a": [0.10, 0.30, 0.60], "b": [0.60, 0.10, 0.30]},
        {"a": [0.05, 0.15, 0.80], "b": [0.15, 0.80, 0.05]},
    ]

    axiom2_samples = [
        {"a": [0.25, 0.75], "b": [0.25, 0.75, 0.00]},
        {"a": [0.40, 0.60], "b": [0.40, 0.60, 0.00]},
        {"a": [0.10, 0.30, 0.60], "b": [0.10, 0.30, 0.60, 0.00]},
        {"a": [0.05, 0.15, 0.80], "b": [0.05, 0.15, 0.80, 0.00]},
    ]

    axiom3_samples = [
        {"a": [0.25, 0.75], "b": [0.30, 0.70]},
        {"a": [0.40, 0.60], "b": [0.41, 0.59]},
        {"a": [0.10, 0.30, 0.60], "b": [0.12, 0.30, 0.58]},
        {"a": [0.05, 0.15, 0.80], "b": [0.10, 0.15, 0.75]},
    ]

    axiom4_samples = []
    for i in range(2, 6):
        new_list = [1.0/i for j in range(i)]
        axiom4_samples.append(new_list)

    if __debug__:
        print("axiom4_samples:", axiom4_samples)

    axiom5_samples = [
        [0.25, 0.75],
        [0.40, 0.60],
        [0.10, 0.30, 0.60],
        [0.05, 0.15, 0.80],
    ]

    if CONVERT_PROPORTION_TO_ABUNDANCE:
        print(f"Converting to abundance with ratio {CONVERT_PROPORTION_TO_ABUNDANCE_RATIO}")
        axiom1_samples = [{"a": to_abundance(i["a"]), "b": to_abundance(i["b"])} for i in axiom1_samples]
        axiom2_samples = [{"a": to_abundance(i["a"]), "b": to_abundance(i["b"])} for i in axiom2_samples]
        axiom3_samples = [{"a": to_abundance(i["a"]), "b": to_abundance(i["b"])} for i in axiom3_samples]
        axiom4_samples = [to_abundance(i) for i in axiom4_samples]
        axiom5_samples = [to_abundance(i) for i in axiom5_samples]
    
    d2_samples = copy.deepcopy(axiom5_samples)
    d3_samples = copy.deepcopy(axiom5_samples)
    d4_samples = copy.deepcopy(axiom5_samples)
    p1_samples = copy.deepcopy(axiom1_samples)
    for i in range(len(p1_samples)):
        p1_samples[i]["b"][0] = math.inf
    p2_samples = copy.deepcopy(axiom2_samples)

    sw2_samples = [
        {"a": [80, 40, 20, 10, 1], "b": [80, 40, 20, 10, 0.5]} # from original paper, 2a and 2b
    ]
    sw3_samples = [
        {"a": [80, 40, 20, 10], "b": [80, 40, 20, 10, 0.5]} # from original paper, 3a and 3b
    ]
    # is sw5 testable automatically?
    sw6_samples = [
        [375, 375, 375, 375] # from original paper, 6
    ]
    sw7_samples = [ # from original paper, 7a, 7b, 7c, 7d, 7e, 7f
        {
            "a": [999, 1],
            "b": [900, 100],
            "c": [800, 200],
            "d": [700, 300],
            "e": [600, 400],
            "f": [500, 500],
        }
    ]
    # is sw8 testable automatically?
    # is sw9 testable automatically?
    # is sw10 testable automatically?
    sw11_samples = [
        {"a": [600, 450, 300, 150], "b": [800, 400, 200, 100]} # from original paper, 11a and 11b
    ] # is sw11 testable automatically?
    # is sw12 testable automatically?
    sw13_samples = [
        {"a": [1000, 1000, 1000, 1], "b": [1000, 1, 1, 1]} # from original paper, 13a and 13b
    ]
    sw14_samples = [
            {"a": [1000, 1, 1, 1, 1, 1], "b": [1000, 1000, 1000, 1, 1, 1], "c": [1000, 1000, 1000, 1000, 1000, 1]} # from original paper, 14a 14b and 14c
    ]

    functions = []
    axiom1_result_list = []
    axiom2_result_list = []
    axiom3_result_list = []
    axiom4_result_list = []
    axiom5_result_list = []
    d1_result_list = []
    d2_result_list = []
    d3_result_list = []
    d4_result_list = []
    p1_result_list = []
    p2_result_list = []
    sw2_result_list = []
    sw3_result_list = []
    sw6_result_list = []
    sw7_result_list = []
    sw13_result_list = []
    sw14_result_list = []

    for i in sorted(list(set(globals()) | set(locals()))):
        if not i.startswith("df_"):
            continue
        print(i.center(64, "="))

        functions.append(i)

        axiom1_result = 1
        for j in axiom1_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            if ROUND_TO_CORRECT_FP_INACCURACY_AXIOM1 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if local_result_a != local_result_b:
                axiom1_result = 0
                print(f"axiom1: {i} failed on {j}; {local_result_a} != {local_result_b}")
                break
        axiom1_result_list.append(axiom1_result)

        axiom2_result = 1
        for j in axiom2_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            if ROUND_TO_CORRECT_FP_INACCURACY_AXIOM2 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if local_result_a != local_result_b:
                axiom2_result = 0
                print(f"axiom2: {i} failed on {j}; {local_result_a} != {local_result_b}")
                break
        axiom2_result_list.append(axiom2_result)

        axiom3_result = 1
        for j in axiom3_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            if ROUND_TO_CORRECT_FP_INACCURACY_AXIOM3 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if local_result_a > local_result_b:
                axiom3_result = 0
                print(f"axiom3: {i} failed on {j}; {local_result_a} > {local_result_b}")
                break
        axiom3_result_list.append(axiom3_result)

        axiom4_result = 1
        for j in axiom4_samples:
            local_result_a = globals()[i](j[:])
            if ROUND_TO_CORRECT_FP_INACCURACY_AXIOM4 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            #print(f"axiom4 log: {local_result_a} != {len(j)} (??)")
            if local_result_a != len(j):
                axiom4_result = 0
                print(f"axiom4: {i} failed on {j}; {local_result_a} != {len(j)}")
                break
        axiom4_result_list.append(axiom4_result)

        axiom5_result = 1
        breaking = False
        for j in axiom5_samples:
            for k in range(2,5):
                local_result_a = globals()[i](j[:])
                j_alt = j * k
                for m in range(len(j_alt)):
                    j_alt[m] = j_alt[m] / k
                local_result_b = globals()[i](j_alt[:])
                local_result_a_times_k = local_result_a * k
                if ROUND_TO_CORRECT_FP_INACCURACY_AXIOM5 == True:
                    local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                    local_result_a_times_k = round(local_result_a_times_k, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                    local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                if local_result_a_times_k != local_result_b:
                    axiom5_result = 0
                    breaking = True
                    print(f"axiom5: {i} failed on {j}; ({local_result_a}) * {k} != {local_result_b}")
                    break
            if breaking:
                break
        axiom5_result_list.append(axiom5_result)

        d2_result = 1
        breaking = False
        for j in d2_samples:
            for k in range(200,500,25):
                alpha = k / 100
                local_result_a = globals()[i](j[:])
                j_alt = [m * alpha for m in j]
                #removed local ratio from axiom5 because not necessary
                local_result_b = globals()[i](j_alt[:])
                if ROUND_TO_CORRECT_FP_INACCURACY_D2 == True:
                    local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                    local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                if local_result_a != local_result_b:
                    d2_result = 0
                    breaking = True
                    print(f"d2: {i} failed on {j}; {local_result_a} != {local_result_b}")
                    break
            if breaking:
                break
        d2_result_list.append(d2_result)

        d3_result = 1
        breaking = False
        for j in d3_samples:
            for k in range(200,500,25):
                alpha = k / 100
                local_result_a = globals()[i](j[:])
                j_alt = [m + alpha for m in j]
                #removed local ratio from axiom5 because not necessary if ALLOW_RATIO == True
                local_result_b = globals()[i](j_alt[:])
                if ROUND_TO_CORRECT_FP_INACCURACY_D3 == True:
                    local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                    local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                if local_result_a >= local_result_b:
                    d3_result = 0
                    breaking = True
                    print(f"d3: {i} failed on {j}; {local_result_a} != {local_result_b}")
                    break
            if breaking:
                break
        d3_result_list.append(d3_result)

        d4_result = 1
        breaking = False
        for j in d4_samples:
            for k in range(2,5):
                local_result_a = globals()[i](j[:])
                j_alt = j * k
                for m in range(len(j_alt)):
                    j_alt[m] = j_alt[m] / k
                local_result_b = globals()[i](j_alt[:])
                if ROUND_TO_CORRECT_FP_INACCURACY_D4 == True:
                    local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                    local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                if local_result_a != local_result_b:
                    d4_result = 0
                    breaking = True
                    print(f"d4: {i} failed on {j}; {local_result_a} != {local_result_b}")
                    break
            if breaking:
                break
        d4_result_list.append(d4_result)


        p1_result = 1
        for j in p1_samples:
            local_result_a = globals()[i](j["a"])
            try:
                local_result_b = globals()[i](j["b"])
            except Exception as e:
                print(e)
                print(f"Skipping P1 testing for {i} because of infinity-related issue")
                p1_result = 0
                break
            if ROUND_TO_CORRECT_FP_INACCURACY_P1 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if not (local_result_a > 0 and local_result_b == 0): #is it really correct to say that as one approaches inf, equity reaches 0?
                p1_result = 0
                print(f"p1: {i} failed on {j}; not ({local_result_a} > 0 and {local_result_b} == 0)")
                break
        p1_result_list.append(p1_result)


        p2_result = 1
        for j in p2_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            if ROUND_TO_CORRECT_FP_INACCURACY_P2 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            # is this correct? careful because of inverse meaning
            if not (local_result_b < local_result_a):
                p2_result = 0
                print(f"p2: {i} failed on {j}; {local_result_a} <= {local_result_b}")
                break
        p2_result_list.append(p2_result)

        sw2_result = 1
        for j in sw2_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            if ROUND_TO_CORRECT_FP_INACCURACY_SW2 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if not (local_result_a > local_result_b):
                sw2_result = 0
                print(f"sw2: {i} failed on {j}; {local_result_a} <= {local_result_b}")
                break
        sw2_result_list.append(sw2_result)

        sw3_result = 1
        for j in sw3_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            if ROUND_TO_CORRECT_FP_INACCURACY_SW3 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if not (local_result_a > local_result_b):
                sw3_result = 0
                print(f"sw3: {i} failed on {j}; {local_result_a} <= {local_result_b}")
                break
        sw3_result_list.append(sw3_result)

        sw6_result = 1
        for j in sw6_samples:
            local_result_a = globals()[i](j)
            if ROUND_TO_CORRECT_FP_INACCURACY_SW6 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if not (local_result_a == 1.0):
                sw6_result = 0
                print(f"sw6: {i} failed on {j}; {local_result_a} != 1.0")
                break
        sw6_result_list.append(sw6_result)

        sw7_result = 1
        for j in sw7_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            local_result_c = globals()[i](j["c"])
            local_result_d = globals()[i](j["d"])
            local_result_e = globals()[i](j["e"])
            local_result_f = globals()[i](j["f"])
            if ROUND_TO_CORRECT_FP_INACCURACY_SW7 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_c = round(local_result_c, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_d = round(local_result_d, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_e = round(local_result_e, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_f = round(local_result_f, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if not (local_result_a < local_result_b and local_result_a < local_result_c and local_result_a < local_result_d and local_result_a < local_result_e and local_result_a < local_result_f):
                sw7_result = 0
                print(f"sw7: {i} failed on {j}; {local_result_a} is not the minimum value: {[local_result_a, local_result_b, local_result_c, local_result_d, local_result_e, local_result_f]}")
                break
        sw7_result_list.append(sw7_result)

        sw13_result = 1
        for j in sw13_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            if ROUND_TO_CORRECT_FP_INACCURACY_SW13 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if not (local_result_a == local_result_b):
                sw13_result = 0
                print(f"sw13: {i} failed on {j}; {local_result_a} != {local_result_b}")
                break
        sw13_result_list.append(sw13_result)

        sw14_result = 1
        for j in sw14_samples:
            local_result_a = globals()[i](j["a"])
            local_result_b = globals()[i](j["b"])
            local_result_c = globals()[i](j["c"])
            if ROUND_TO_CORRECT_FP_INACCURACY_SW14 == True:
                local_result_a = round(local_result_a, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_b = round(local_result_b, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
                local_result_c = round(local_result_c, ROUND_TO_CORRECT_FP_INACCURACY_NDIGITS)
            if not (local_result_b > local_result_a and local_result_b > local_result_c):
                sw14_result = 0
                print(f"sw14: {i} failed on {j}; {local_result_b} <= {local_result_a} or {local_result_b} <= {local_result_c}")
                break
        sw14_result_list.append(sw14_result)

    df = pd.DataFrame(
        {
            "function": functions,
            "a1": axiom1_result_list,
            "a2": axiom2_result_list,
            "a3/d1": axiom3_result_list,
            "a4": axiom4_result_list,
            "a5": axiom5_result_list,
            "d2/sw4": d2_result_list,
            "d3": d3_result_list,
            "d4/sw1": d4_result_list,
            "p1/sw7~": p1_result_list,
            "p2": p2_result_list,
            "sw2": sw2_result_list,
            "sw3": sw3_result_list,
            "sw6": sw6_result_list,
            "sw7": sw7_result_list,
            "sw13": sw13_result_list,
            "sw14": sw14_result_list,
        }
    )
    print("Morales et al. (2020):")
    print(f"{'a1':<8} Symmetry")
    print(f"{'a2':<8} Expansibility")
    print(f"{'a3':<8} Transfer Principle")
    print(f"{'a4':<8} Normalization")
    print(f"{'a5':<8} Replication")
    print("Hurley and Rickard (2009):")
    print(f"{'d1':<8} Robin Hood")
    print(f"{'d2':<8} Scaling")
    print(f"{'d3':<8} Rising Tide")
    print(f"{'d4':<8} Cloning")
    print(f"{'p1':<8} Bill Gates")
    print(f"{'p2':<8} Babies")
    print("Smith and Wilson (1996):")
    print(f"{'sw1':<8} <-> d4")
    print(f"{'sw2':<8} Slight increase in lowest type increases diversity")
    print(f"{'sw3':<8} Adding new lowest type decreases diversity")
    print(f"{'sw4':<8} <-> d2")
    print(f"{'sw6':<8} Maximum value (i.e., when perfectly balanced) equals 1")
    print(f"{'sw7':<8} Minimal value when as unequal as it gets")
    print(f"{'sw13':<8} \"Symmetric with regard to minor and abundant species\"")
    print(f"{'sw14':<8} \"Skewed distributions should give a lower value\"")

    n_rows = len(df)
    n_cols = len(df.columns) - 1
    print(n_rows, n_cols)

    df.index = df['function']
    df.pop("function")
    df['ratio'] = df.sum(axis=1, numeric_only=True)
    df['ratio'] /= n_cols
    df_lower_sum = pd.DataFrame([df.sum(axis=0, numeric_only=True)])
    df_lower_sum /= n_rows
    df_lower_sum.index = ["ratio"]
    df = pd.concat((df,df_lower_sum))
    
    #df['n_axioms'] = df['a1'] + df['a2'] + df['a3/d1'] + df['a4'] + df['a5'] + df['d2/sw4'] + df['d3'] + df['d4/sw1'] + df['p1/sw7~'] + df['p2'] + df['sw2'] + df['sw3'] + df['sw6'] + df['sw7'] + df['sw13'] + df['sw14'] # DO NOT REMOVE

    print(df)
    df.to_csv("output_axioms_df.csv", sep=";", encoding="utf-8")

    #### INTRASET ####

    REFERENCE_TOKENS = ['great', 'nice', 'good', 'bad']
    functions_ = []
    result_list = {
        "s1": [],
        "s2": [],
        "s3": [],
        "s4": [],
        "s5": [],
    }

    for i in sorted(list(set(globals()) | set(locals()))):
        if not "intraset_disparity" in i:
            continue

        functions_.append(i)

        # s1
        tokens = REFERENCE_TOKENS[:1]
        amounts = [5]
        f_result = globals()[i](amounts, tokens)
        result_list["s1"].append(int(f_result == 0.0))

        # s2
        previous_result = None
        success = 1
        for j in range(1, 10):
            tokens = REFERENCE_TOKENS[:1] * j # should it be [:2] or [:1]?
            amounts = [10 for k in range(len(tokens))]
            f_result = globals()[i](amounts, tokens)
            if previous_result != None and f_result < previous_result:
                success = 0
                break
            previous_result = f_result
        result_list["s2"].append(success)

        # s3
        previous_result = None
        success = 1
        n = 10
        for j in range(0, n):
            tokens = REFERENCE_TOKENS[:2] # should it be [:1] instead?
            amounts = [n-1, j]
            f_result = globals()[i](amounts, tokens)
            if previous_result != None and f_result < previous_result:
                success = 0
                break
            previous_result = f_result
        result_list["s3"].append(success)

        # s4
        previous_result = None
        success = 1
        step = 5
        for j in range(step, step*5, step):
            tokens = REFERENCE_TOKENS[:2]
            amounts = [math.floor((1.0/step)*j), math.floor((1.0-(1.0/step))*j)]
            f_result = globals()[i](amounts, tokens, vector_ratio=(j/step))
            if previous_result != None and f_result < previous_result:
                success = 0
                break
            previous_result = f_result
        result_list["s4"].append(success)

        # s5
        success = 1
        n = 5
        for j in range(len(REFERENCE_TOKENS)):
            tokens = [REFERENCE_TOKENS[j]] * n
            amounts = [10 for k in range(len(tokens))]
            f_result = globals()[i](amounts, tokens)
            if f_result != 0.0:
                success = 0
                break
        result_list["s5"].append(success)

    df = pd.DataFrame(result_list)
    df['n_axioms'] = [0] * len(functions_)
    for k in result_list:
        df['n_axioms'] += df[k]
    df.index = functions_

    print("Stirling (2007):")
    print(f"{'s1':<8} \"Scaling of variety\"")
    print(f"{'s2':<8} \"Monotonicity of variety\"")
    print(f"{'s3':<8} \"Monotonicity of balance\"")
    print(f"{'s4':<8} \"Monotonicity of disparity\"")
    print(f"{'s5':<8} \"Scaling of disparity\"")

    print(df)

    df.to_csv("output_axioms_disparity.csv", sep=";", encoding="utf-8")

    return 0

if __name__ == "__main__":
    exit(main())


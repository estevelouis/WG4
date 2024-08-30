import diversutils
from typing import *

def _wrapper_df(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, measure_key: str = "", alpha: float = 1.0, beta: float = 1.0) -> Tuple[float]:
    if (not measure_key.startswith("DF_")) or measure_key not in dir(diversutils):
        raise Exception(f"Undefined diversity function: {measure_key}")

    IS_DISPARITY = measure_key.startswith("DF_DISPARITY")


    if IS_DISPARITY:
        if w2v_path == "":
            raise Exception(f"A w2v_path must be provided.")
        if not os.path.exists(w2v_path):
            raise Exception(f"A w2v_path was given ('{w2v_path}') but it does not exist.")

        graph_index = diversutils.create_empty_graph(0, w2v_num_dimensions)
        if graph_index = -1:
            raise Exception(f"Failed to create graph")
        w2v_index = diversutils.load_w2v(w2v_path)
        if w2v_index = -1:
            raise Exception(f"Failed to load w2v")
        if diversutils.bind_w2v(graph_index, w2v_index):
            raise Exception(f"Failed to bind the vector space to the graph")
    else:
        graph_index = diversutils.create_empty_graph(0, 0)
        if graph_index = -1:
            raise Exception(f"Failed to create graph")


    for key, value in category_dict.items():
        if type(key) != str:
            raise Exception(f"Key of category_dict must of type str. Received type: {type(key)}")
        if type(value) != int:
            value = int(value)

        if value <= 0:
            raise Exception(f"Each category must have a least one element.")

        if IS_DISPARITY:
            add_node_res = diversutils.add_node(graph_index, value, key)
        else:
            add_node_res = diversutils.add_node(graph_index, value)
        if add_node_res != 0:
            raise Exception(f"Failed to add node: {key}, {value}")

    diversutils.compute_relative_proportion(graph_index)

    results = diversutils.individual_measure(graph_index, getattr(diversutils, measure_key), alpha, beta)

    diversutils.free_graph(graph_index)

    if IS_DISPARITY:
        diversutils.free_w2v(w2v_index)

    return results


def DISPARITY_CHAO_ET_AL_FUNCTIONAL(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_DISPARITY_CHAO_ET_AL_FUNCTIONAL')
def DISPARITY_LEINSTER_COBBOLD(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_DISPARITY_LEINSTER_COBBOLD')
def DISPARITY_PAIRWISE(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_DISPARITY_PAIRWISE')
def DISPARITY_RICOTTA_SZEIDL(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_DISPARITY_RICOTTA_SZEIDL')
def DISPARITY_SCHEINER(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_DISPARITY_SCHEINER')
def DISPARITY_STIRLING(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_DISPARITY_STIRLING')
def ENTROPY_GOOD(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_ENTROPY_GOOD')
def ENTROPY_PATIL_TAILLIE(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_ENTROPY_PATIL_TAILLIE')
def ENTROPY_Q_LOGARITHMIC(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_ENTROPY_Q_LOGARITHMIC')
def ENTROPY_RENYI(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_ENTROPY_RENYI')
def ENTROPY_SHANNON_WEAVER(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_ENTROPY_SHANNON_WEAVER')
def INDEX_BERGER_PARKER(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_BERGER_PARKER')
def INDEX_BRILLOUIN(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_BRILLOUIN')
def INDEX_E_BULLA1994(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_E_BULLA1994')
def INDEX_E_HEIP(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_E_HEIP')
def INDEX_E_MCI_PIELOU1969(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_E_MCI_PIELOU1969')
def INDEX_E_MINUS_LN_D_PIELOU1977(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_E_MINUS_LN_D_PIELOU1977')
def INDEX_E_PRIME_CAMARGO1993(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_E_PRIME_CAMARGO1993')
def INDEX_E_VAR_SMITH_AND_WILSON1996(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_E_VAR_SMITH_AND_WILSON1996')
def INDEX_F_2_1_ALATALO1981(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_F_2_1_ALATALO1981')
def INDEX_G_2_1_MOLINARI1989(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_G_2_1_MOLINARI1989')
def INDEX_HILL_EVENNESS(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_HILL_EVENNESS')
def INDEX_JUNGE1994_PAGE22(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_JUNGE1994_PAGE22')
def INDEX_MCINTOSH(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_MCINTOSH')
def INDEX_ONE_MINUS_D(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_ONE_MINUS_D')
def INDEX_ONE_OVER_D_WILLIAMS1964(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_ONE_OVER_D_WILLIAMS1964')
def INDEX_O_BULLA1994(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_O_BULLA1994')
def INDEX_RICHNESS(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_RICHNESS')
def INDEX_SHANNON_EVENNESS(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_SHANNON_EVENNESS')
def INDEX_SIMPSON(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_SIMPSON')
def INDEX_SIMPSON_DOMINANCE(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_SIMPSON_DOMINANCE')
def INDEX_SPECIES_COUNT(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_SPECIES_COUNT')
def INDEX_TYPE_TOKEN_RATIO(category_dict: Dict[str, int], w2v_path: str = "", w2v_num_dimensions: int = 100, alpha: float = 1.0, beta: float = 1.0) -> float:
    return _wrapper_df(category_dict=category_dict, w2v_path=w2v_path, w2v_num_dimensions=w2v_num_dimensions, alpha=alpha, beta=beta, measure_key='DF_INDEX_TYPE_TOKEN_RATIO')



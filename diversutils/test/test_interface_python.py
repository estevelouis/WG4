import diversutils
import math

def test_import_diversutils():
    assert ("diversutils" in globals()), "failed to import diversutils"

def test_create_graph():
    g_index = diversutils.create_empty_graph(0, 0)
    assert (diversutils.free_graph(g_index) == 0), "failed to free graph"

def test_add_node_no_key():
    g_index = diversutils.create_empty_graph(0, 0)
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.free_graph(g_index) == 0), "failed to free graph"

def test_compute_relative_frequencies():
    g_index = diversutils.create_empty_graph(0, 0)
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.compute_relative_proportion(g_index) == 0), "failed to compute relative frequencies"
    assert (diversutils.free_graph(g_index) == 0), "failed to free graph"

def test_shannon_weaver_entropy():
    g_index = diversutils.create_empty_graph(0, 0)
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.compute_relative_proportion(g_index) == 0), "failed to compute relative frequencies"
    assert (diversutils.individual_measure(g_index, diversutils.DF_ENTROPY_SHANNON_WEAVER) == (2.0 * (0.5 * -math.log(0.5, math.e)), math.e ** (2.0 * (0.5 * -math.log(0.5, math.e))))), "failed to have correct result for Shannon-Weaver entropy"
    assert (diversutils.free_graph(g_index) == 0), "failed to free graph"

def test_indices_not_nan():
    g_index = diversutils.create_empty_graph(0, 0)
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.add_node(g_index, 4) == 0), "failed to add node"
    assert (diversutils.compute_relative_proportion(g_index) == 0), "failed to compute relative frequencies"
    df = list(filter(lambda x: x.startswith(("DF_ENTROPY", "DF_INDEX")), dir(diversutils)))
    for i in df:
        res = diversutils.individual_measure(g_index, getattr(diversutils, i))
        print(i, res)
        assert (math.nan not in res), f"nan found in result of {i}"
    assert (diversutils.free_graph(g_index) == 0), "failed to free graph"

"""
def test_disparity_not_nan():
    #g_index = diversutils.create_empty_graph(0, 200)
    g_index = diversutils.create_empty_graph(0, 100)
    #w2v_index = diversutils.load_w2v("data/vectors/vec.bin")
    #w2v_index = diversutils.load_w2v("data/vectors/vec2.bin")
    w2v_index = diversutils.load_w2v("/home/esteve/Documents/thesis/other_repos/word2vec/bin/MWE_S2S_IT_11GB_100d_skip-gram.bin")
    diversutils.bind_w2v(g_index, w2v_index)
    assert (diversutils.add_node(g_index, 4, "vero") == 0), "failed to add node"
    assert (diversutils.add_node(g_index, 4, "altro") == 0), "failed to add node"
    #assert (diversutils.add_node(g_index, 4, "good") == 0), "failed to add node"
    #assert (diversutils.add_node(g_index, 4, "bad") == 0), "failed to add node"
    #assert (diversutils.add_node(g_index, 4, "eight") == 0), "failed to add node"
    #assert (diversutils.add_node(g_index, 4, "from") == 0), "failed to add node"
    assert (diversutils.compute_relative_proportion(g_index) == 0), "failed to compute relative frequencies"
    df = list(filter(lambda x: x.startswith(("DF_DISPARITY",)), dir(diversutils)))
    for i in df:
        res = diversutils.individual_measure(g_index, getattr(diversutils, i))
        print(i, res)
        assert (math.nan not in res), f"nan found in result of {i}"
    assert (diversutils.free_graph(g_index) == 0), "failed to free graph"
    assert (diversutils.free_w2v(w2v_index) == 0), "failed to free w2v"
"""

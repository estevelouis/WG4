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


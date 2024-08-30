from _diversutils import _attach_distance_matrix
from _diversutils import *
import numpy as _np
import ctypes as _ct

def attach_distance_matrix(graph_index: int, matrix: _np.ndarray, fp_mode: int = FP32):
    if type(graph_index) != int:
        raise Exception("graph_index must be of type int")
    if type(matrix) != _np.ndarray:
        raise Exception("matrix must be of type numpy.ndarray")
    if type(fp_mode) != int:
        raise Exception("fp_mode must be of type int")
    if fp_mode == FP32:
        matrix = matrix.astype(_np.float32)
    elif fp_mode == FP64:
        matrix = matrix.astype(_np.float64)
    else:
        raise Exception("fp_mode must be either FP32 or FP64")

    bfr = matrix.tobytes()
    #void_p = _ct.addressof(bfr)
    c_bfr = _ct.c_buffer(bfr)
    void_p = _ct.addressof(c_bfr)

    _attach_distance_matrix(graph_index, void_p, fp_mode)


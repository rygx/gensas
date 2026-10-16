import numpy as np

from gensas_server.results import reduce_indices


def test_small_results_are_kept_whole():
    t = np.linspace(0, 1, 10)
    assert list(reduce_indices(t, [t], 10)) == list(range(10))
    assert list(reduce_indices(t, [t], None)) == list(range(10))


def test_reduction_keeps_ends_events_and_extremes():
    t = np.linspace(0, 10, 2001)
    t[1000:] = np.concatenate(([t[1000]], t[1000:-1]))  # two samples at one time: an event pair
    t[1000] = t[1001]
    x = np.sin(t)
    x[777] = 5.0  # a narrow spike
    idx = reduce_indices(t, [x], 100)
    assert idx[0] == 0 and idx[-1] == 2000
    assert 777 in idx
    assert 1000 in idx and 1001 in idx
    assert len(idx) <= 130
    assert list(idx) == sorted(set(idx))


def test_nan_columns_do_not_break_reduction():
    t = np.linspace(0, 1, 500)
    x = np.full(500, np.nan)
    idx = reduce_indices(t, [x], 50)
    assert len(idx) >= 2

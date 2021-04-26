from bgmol.systems import ChignolinC22Implicit


def test_chignolin_tics():
    cgn = ChignolinC22Implicit()
    tics = cgn.to_tics(cgn.positions[None, ...], eigs_kept=10)
    assert tics.shape == (1, 10)
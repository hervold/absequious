from absequious import utils


def test_translate_six():
    assert utils.translate_six("ATGCAACCGATGACAAA") == [
        (utils.Dir.fwd, 0, "MQPMT"),
        (utils.Dir.fwd, 1, "CNR*Q"),
        (utils.Dir.fwd, 2, "ATDDK"),
        (utils.Dir.rev_comp, 0, "FVIGC"),
        (utils.Dir.rev_comp, 1, "LSSVA"),
        (utils.Dir.rev_comp, 2, "CHRLH"),
    ]

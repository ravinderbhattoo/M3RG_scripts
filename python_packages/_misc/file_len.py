def len(fname):
    r"""
    Get length of file.
    """

    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

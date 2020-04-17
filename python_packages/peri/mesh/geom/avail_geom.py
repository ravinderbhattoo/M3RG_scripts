r"""
    This is a misc module.
    Please see the content and functions.
    Keep adding new ones :)
    """

def list_all():
    from os.path import dirname, basename, isfile
    import glob
    modules = glob.glob(dirname(__file__)+"/*.py")
    list_l1 = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
    list_l1.remove('avail_geom')
    return list_l1
    import setup_halos as sh
    import get_stats as gs
    import numpy as np

    halocat, HODmodel = sh.setup_halos()
    bins = np.logspace(np.log10(5.0), np.log10(15.0), 15)
    bcens, wtheta = gs.get_stats(halocat, HODmodel, bins)

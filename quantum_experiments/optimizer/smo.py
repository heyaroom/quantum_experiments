from __future__ import division
import numpy as np
import copy as copy
from scipy.optimize import minimize, OptimizeResult

def optimize(model, p_seed, iteration, initp=None):
    n_param = model.n_param

    np.random.seed(p_seed)
    if initp is None:
        param = np.random.random(size=n_param) * 2 * np.pi
    else:
        param = initp

    maxfev = 1 + 2*n_param*iteration

    res  = minimize(
        model.step,
        copy.copy(param),
        options={'maxfev': maxfev,"reset_interval":-1},
        method=nakanishi_fujii_todo,
        callback=model.callback
    )

def nakanishi_fujii_todo(fun, x0, args=(), maxfev=1024, reset_interval=32, eps=1e-32, callback=None, **_):
    x0 = np.asarray(x0)
    recycle_z0 = None
    niter = 0
    funcalls = 0

    while True:

        idx = niter % x0.size

        if reset_interval > 0:
            if niter % reset_interval == 0:
                recycle_z0 = None

        if recycle_z0 is None:
            z0 = fun(np.copy(x0), *args)
            funcalls += 1
        else:
            z0 = recycle_z0

        p = np.copy(x0)
        p[idx] = x0[idx] + np.pi / 2
        z1 = fun(p, *args)
        funcalls += 1

        p = np.copy(x0)
        p[idx] = x0[idx] - np.pi / 2
        z3 = fun(p, *args)
        funcalls += 1

        z2 = z1 + z3 - z0
        c = (z1 + z3) / 2
        a = np.sqrt((z0 - z2) ** 2 + (z1 - z3) ** 2) / 2
        b = np.arctan((z1 - z3) / ((z0 - z2) + 1e-32 * (z0 == z2))) + x0[idx]
        b += 0.5 * np.pi + 0.5 * np.pi * np.sign((z0 - z2) + eps * (z0 == z2))

        x0[idx] = b
        recycle_z0 = c - a

        if callback is not None:
            callback(np.copy(x0), np.copy(recycle_z0))

        if funcalls >= maxfev:
            break

        niter += 1

    return OptimizeResult(fun=fun(np.copy(x0)), x=x0, nit=niter, nfev=funcalls, success=(niter > 1))
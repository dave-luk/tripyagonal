import sympy as sp
from tripyagonal import AbstractMatrix


def latex(obj):
    try:
        if isinstance(obj, list):
            l = []
            for e in obj:
                l.append(latex(e))
            return l
        elif isinstance(obj, AbstractMatrix):
            return obj.latex()
        return sp.latex(obj)
    except:
        raise


def evalf(obj):
    try:
        if isinstance(obj, list):
            l = []
            for e in obj:
                l.append(evalf(e))
            return l
        elif isinstance(obj, AbstractMatrix):
            return obj.evalf()
        return obj.evalf()
    except:
        raise


def pretty(obj):
    try:
        if isinstance(obj, list):
            l = []
            for e in obj:
                l.append(pretty(e))
            return l
        elif isinstance(obj, AbstractMatrix):
            return obj.pretty()
        return sp.pretty(obj)
    except:
        raise


def pprint(obj):
    print(pretty(obj))
    return

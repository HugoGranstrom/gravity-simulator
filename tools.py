import platform

### containerVector ###
from collections import namedtuple

conVec = namedtuple("conVec", "x y")

### Vector Class ###
try:
    if platform.python_implementation() == "PyPy":
        import vector  # use pure python vector for PyPy

        vector = vector.vector
    else:
        import cyvector.vector

        vector = cyvector.vector
except:
    import cyvector

    vector = cyvector.vector

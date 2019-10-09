import faulthandler
faulthandler.enable()
import numpy as np

if True: # "scope"
    print("Rx, no settings")
    import qop
    c = qop.Circuit(1)
    x = qop.Gate('rx', parameters=[0.1])
    c.add_gate(x, 0)
    result = c.optimize([[1, 0], [0, -1]])
    print(result)

if True:
    print("Rx, with settings")
    import qop
    c = qop.Circuit(1)
    x = qop.Gate('rx', parameters=[0.1])
    c.add_gate(x, 0)
    result = c.optimize([[1, 0], [0, -1]], settings={
        'ada': {
            'rho': 0.9,
            'epsilon': 0.001,
        }
    })
    print(result)
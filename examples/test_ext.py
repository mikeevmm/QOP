import numpy as np
import faulthandler
faulthandler.enable()
""" 
print("Rx, no settings")
import qop
c = qop.Circuit(1)
x = qop.Gate('rx', parameters=[0.1])
c.add_gate(x, 0)
result = c.optimize([[1, 0], [0, -1]])
print(result)

print("Rx, with settings")
import qop
c = qop.Circuit(1)
x = qop.Gate('rx', parameters=[0.1])
c.add_gate(x, 0)
result = c.optimize([[1, 0], [0, -1]], settings={
    'ada': {
        'rho': 0.9,
        'epsilon': 0.001,
    },
    'optimize': {
        'max_iterations': 5
    }
})
print(result)
 """
print("Multiple")
import qop
c = qop.Circuit(2)
x1 = qop.Gate('rx', parameters=[0.1])
x2 = qop.Gate('rx', parameters=[0.1])
c.add_gate(x1, 0)
c.add_gate(x2, 1)
result = c.optimize([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, -1]],
                    settings={
    'ada': {
        'rho': 0.95
    },
    'optimize': {
        'gates': [x1],
        'deltas': [[0.01]],
        'max_iterations': 1e4
    }
})
print(result)
print("What gates?")
print(c.get_gates())
print("x2 is unaffected: ", x2.get_parameters())
print("and we can reparameterize it")
x2.reparameterize([2])
print(x2.get_parameters())
print(x2.get_matrix())

print("Now let's run the circuit once...")
in_state = [1, 0, 0, 0]
out = c.run(in_state)
print(out)
import qop
import math
import random
import numpy as np

qubit_count = 3

circuit = qop.Circuit(qubit_count)
for line in range(qubit_count):
    h = qop.Gate('h')
    circuit.add_gate(h, line)
    for n in range(2, 2 + qubit_count - 1 - line):
        rz = qop.Gate('rz', parameters=[2*math.pi/2**n])
        circuit.add_gate(rz, line, control=line + n - 1)

# Pick a random number
x = random.randint(0, 2**qubit_count - 1)
print(f"in state is {x}")
in_state = [1. if u == x else 0. for u in range(2**qubit_count)]
print("Two norm before:")
print(sum(in_state[i].conjugate()*in_state[i]
          for i in range(2**qubit_count)))

out_state = circuit.run(in_state)
print(out_state)

print("Two norm after:")
print(sum(out_state[i].conjugate()*out_state[i]
          for i in range(2**qubit_count)))

print("Verifying coefficients")
for i, coef in enumerate(out_state):
    predicted = 1/math.sqrt(2**qubit_count) * np.exp(2*math.pi*1.0j*x*i/2**qubit_count)
    if abs(coef - predicted) > 1e-3:
        print(f"Inconsistency of coefficent #{i}")
        print(f"Expected {predicted} got {coef}")
    else:
        print(f"Coefficient {i} ✔")
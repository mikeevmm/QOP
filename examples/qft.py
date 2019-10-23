import qop
import math
import random
import numpy as np

qubit_count = 4

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
bit_mirror = lambda n: int(f'{{:0{qubit_count}b}}'.format(n)[::-1], 2)
got = ((out_state[i]) for i in range(2**qubit_count))

for basis, coef in enumerate(got):
    expected = 1/np.sqrt(2**qubit_count) * np.exp(2*np.pi*1.0j*bit_mirror(x)*basis/2**qubit_count)
    if abs(expected - coef) > 1e-4:
        print("DISCREPANCY")
        print(f"In |{basis}> expected {expected} got {coef}")
    else:
        print("✔")

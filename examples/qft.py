import qop
import math
import random

qubit_count = 5

circuit = qop.Circuit(qubit_count)
for line in range(qubit_count):
    h = qop.Gate('h')
    circuit.add_gate(h, line)
    for n in range(2, 2 + qubit_count - line - 1):
        rz = qop.Gate('rz', parameters=[2*math.pi/2**n])
        circuit.add_gate(rz, line, control=line + n - 1)

# Pick a random number
x = random.randint(0, 2**qubit_count)
in_state = [1. if u == x else 0. for u in range(2**qubit_count)]
print("Two norm:")
print(sum(in_state[i].conjugate()*in_state[i]
          for i in range(2**qubit_count)))

out_state = circuit.run(in_state)
print(out_state)

print("Two norm:")
print(sum(out_state[i].conjugate()*out_state[i]
          for i in range(2**qubit_count)))

import numpy.linalg
import numpy as np
from time import time
import qop
import random

random.seed(8081902)

layer_count = 4
qubit_count = 3
matrix = [[1.50000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,
           0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, ],
          [0.000000000000000, -0.500000000000000,   2.00000000000000,  0.000000000000000,
           2.00000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, ],
          [0.000000000000000,   2.00000000000000, -0.500000000000000,  0.000000000000000,
           2.00000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, ],
          [0.000000000000000,  0.000000000000000,  0.000000000000000, -0.500000000000000,
           0.000000000000000,   2.00000000000000,   2.00000000000000,  0.000000000000000, ],
          [0.000000000000000,   2.00000000000000,   2.00000000000000,  0.000000000000000, -
           0.500000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, ],
          [0.000000000000000,  0.000000000000000,  0.000000000000000,   2.00000000000000,
           0.000000000000000, -0.500000000000000,   2.00000000000000,  0.000000000000000, ],
          [0.000000000000000,  0.000000000000000,  0.000000000000000,   2.00000000000000,
           0.000000000000000,   2.00000000000000, -0.500000000000000,  0.000000000000000, ],
          [0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,   1.50000000000000]]
c = qop.Circuit(qubit_count)
ry_gates = []
for _ in range(layer_count):
    for i in range(0, 1):
        for qubit in range(qubit_count):
            ry = qop.Gate('ry', parameters=[
                          (random.random() - 0.5) * 2 * 3.1415926])
            c.add_gate(ry, qubit)
            ry_gates.append(ry)
        for qubit in range(0, qubit_count, 2):
            z = qop.Gate('z')
            c.add_gate(z, qubit + i, (qubit + i + 1) % qubit_count)
for qubit in range(qubit_count):
    ry = qop.Gate('ry', parameters=[(random.random() - 0.5) * 2 * 3.1415926])
    ry_gates.append(ry)
    c.add_gate(ry, qubit)

start = time()
results, broke = c.optimize(matrix, settings={'ada': {'rho': 0.90}, 'optimize': {
    'stop_at': 1e-4}})
end = time()
print(results)
print(f"in {end - start}s")
print(("" if broke else "not ") + "breaking after max iterations")

print("output:")
out = c.run([1] + [0] * (2**qubit_count - 1))
print(out)

out = np.array(out, dtype=np.complex).T  # Column vector
ham = np.array(matrix, dtype=np.complex)
print("Out state 2-norm")
print(np.linalg.norm(out, ord=2))
print("Hamiltonian contraction")
print(np.conjugate(out.T) @ ham @ out)

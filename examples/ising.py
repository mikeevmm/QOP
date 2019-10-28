import numpy.linalg
import numpy as np
from time import time
import qop
import random

# random.seed(8081902)

layer_count = 4
spin_count = 6

# Build the hamiltonian
sx = np.matrix([[0, 1], [1, 0]], dtype=np.complex)
sy = np.matrix([[0, -1.0j], [1.j, 0]], dtype=np.complex)
sz = np.matrix([[1, 0], [0, -1]], dtype=np.complex)

# Return the matrix for a spin of type u (i.e. x, y, z) at position j,
# for a q qubit spin ring


def suj(q, u, j):
    return np.kron(np.kron(np.eye(2**j), u), np.eye(2**(q - j - 1)))


matrix = sum(suj(spin_count, sx, i) * suj(spin_count, sx, (i + 1) % spin_count) +
             suj(spin_count, sy, i) * suj(spin_count, sy, (i + 1) % spin_count) +
             0.5 * suj(spin_count, sz, i) *
             suj(spin_count, sz, (i + 1) % spin_count)
             for i in range(spin_count))

# Run the simulation
c = qop.Circuit(spin_count)
ry_gates = []
for _ in range(layer_count):
    for i in range(2):
        for qubit in range(spin_count):
            ry = qop.Gate('ry', parameters=[
                (random.random() - 0.5) * 2 * 3.1415926])
            c.add_gate(ry, qubit)
            ry_gates.append(ry)
        for qubit in range(0, spin_count, 2):
            z = qop.Gate('z')
            c.add_gate(z, qubit + i, (qubit + i + 1) % spin_count)
for qubit in range(spin_count):
    ry = qop.Gate('ry', parameters=[(random.random() - 0.5) * 2 * 3.1415926])
    ry_gates.append(ry)
    c.add_gate(ry, qubit)

with open('../test.txt', 'w') as test:
    start = time()
    results, broke = c.optimize(matrix, writer=test)
    out = np.array(c.run([1] + [0] * (2**spin_count - 1)))
    print(np.conjugate(out.T) @ matrix @ out)
    end = time()

print("results:")
print(results)
print(f"in {end - start}s")
print(("" if broke else "not ") + "breaking after max iterations")

print("output:")
out = c.run([1] + [0] * (2**spin_count - 1))

expct_value = 0
for i in range(2**spin_count):
    for j in range(2**spin_count):
        expct_value += out[i].conjugate() * matrix[i, j] * out[j]
print(expct_value)

out = np.array(out, dtype=np.complex).T  # Column vector
print("Out state 2-norm")
print(np.linalg.norm(out, ord=2))
print("Hamiltonian contraction")
print(np.conjugate(out.T) @ matrix @ out)

print("Minimum eigenvalue, from numpy:")
print(min(np.linalg.eigvals(matrix)))

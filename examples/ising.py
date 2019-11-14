import qop
import numpy as np
import random
from functools import reduce

# Parameters
qubit_count = 5
layer_count = 4
settings = {
    'ada': {
        'rho': 0.99,
    },
    'optimize': {
        'max_iterations': 3000
    }
}

# Definitions
sx = np.array([(0, 1),
               (1, 0)])
sz = np.array([(1, 0),
               (0, -1)])

def sui(su, i, N):
    return np.kron(np.kron(reduce(np.kron, [np.eye(2) for _ in range(i)], 1), su), reduce(np.kron, [np.eye(2) for _ in range(i+1,N)], 1))

# Build the hamiltonian
hamiltonian = sum([sui(sx, i, qubit_count) @ sui(sx, (i+1)%qubit_count, qubit_count) for i in range(qubit_count)])

# Build the circuit
c = qop.Circuit(qubit_count)
ry_gates = []
for _ in range(layer_count):
    for i in (0,1):
        for qubit in range(qubit_count):
            ry = qop.Gate('ry', parameters=[
                (random.random() - 0.5) * 2 * 3.1415926])
            c.add_gate(ry, qubit)
            ry_gates.append(ry)
        for qubit in range(0, qubit_count, 2):
            z = qop.Gate('z')
            c.add_gate(z, (qubit + i) % qubit_count, (qubit + i + 1) % qubit_count)
for qubit in range(qubit_count):
    ry = qop.Gate('ry', parameters=[(random.random() - 0.5) * 2 * 3.1415926])
    ry_gates.append(ry)
    c.add_gate(ry, qubit)

# Perform optimization
with open('ising-output.tsv', 'w') as writer:
	results, broke = c.optimize(hamiltonian, settings=settings, writer=writer)

print("output:")
out = c.run([1] + [0] * (2**qubit_count - 1))

expct_value = 0
for i in range(2**qubit_count):
    for j in range(2**qubit_count):
        expct_value += out[i].conjugate() * hamiltonian[i, j] * out[j]
print(expct_value)

out = np.array(out, dtype=np.complex).T  # Column vector
print("Out state 2-norm")
print(np.linalg.norm(out, ord=2))
print("Hamiltonian contraction")
print(np.conjugate(out.T) @ hamiltonian @ out)

print("Minimum eigenvalue, from numpy:")
print(min(np.linalg.eigvals(hamiltonian)))

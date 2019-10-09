import faulthandler
faulthandler.enable()
import numpy as np

import qop
c = qop.Circuit(1)
x = qop.Gate('rx', parameters=[0.1])
c.add_gate(x, 0)
result = c.optimize([[1, 0], [0, -1]])
print(result)
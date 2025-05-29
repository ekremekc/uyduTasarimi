import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
# Arrays to store data
x_vals = []
y_vals = []

# List all files in the current directory
for filename in sorted(os.listdir(".")):
    if filename.startswith("T_diff_") and filename.endswith(".txt"):
        with open(filename) as f:
            lines = f.readlines()
            if len(lines) >= 2:
                x = float(lines[0].strip())
                y = float(lines[1].strip())
                x_vals.append(x)
                y_vals.append(y)

# Plotting
plt.plot(x_vals, y_vals, 'o-', label='T_diff')
plt.xlabel("Thickness (m)")
plt.ylabel("Temperature Difference (C)")
# plt.title("T_diff vs Thickness")
plt.grid(True)
# plt.legend()
plt.gca().invert_xaxis()
plt.tight_layout()
plt.savefig("T_diff.png")
# plt.show()
plt.close()

import numpy as np

def p_contact(T, N_c, f, d_c, A_contact):
    P = T*N_c/(f*d_c*A_contact)
    return P

def calculate_torque(P, N_c, f, d_c, A_contact):
    T = P*(f*d_c*A_contact)/N_c
    return T

# T = np.linspace(1000, 10000, 20) #N
N_c = 4
f = 0.2
d_c = 5 #mm
A_contact = np.pi*d_c**2/4
Ps = np.array([34, 69, 172, 345])*1E-3
Ts = calculate_torque(Ps, N_c, f, d_c, A_contact)
print(Ts)

plt.plot(Ps, Ts, '-s')
plt.xlabel("Contact pressure (N/mm2)")
plt.ylabel("Torque (N)")
plt.tight_layout()

plt.grid()
plt.show()
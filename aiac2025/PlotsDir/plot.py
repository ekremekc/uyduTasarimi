import os
import matplotlib.pyplot as plt

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
plt.xlabel("Thickness")
plt.ylabel("Temperature Difference (C)")
# plt.title("T_diff vs Thickness")
plt.grid(True)
# plt.legend()
plt.gca().invert_xaxis()
plt.tight_layout()
plt.savefig("T_diff.png")
plt.show()
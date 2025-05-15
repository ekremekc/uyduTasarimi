kappa = 167.9 # W/mK 167.9
kappa_pad = 0.2
Q_total = -15  # W
thermal_pad = True


if "__name__" == "__main__":
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_csv("MaterialsDir/al_kappa.txt")

    print(df["Temperature(C)"])
    print(df["ThermalConductivity(microW/mm-dC)"]*1E-3)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(df["Temperature(C)"], df["ThermalConductivity(microW/mm-dC)"]*1E-3, marker='o')
    plt.xlabel("Temperature(C)")
    plt.ylabel("Thermal Conductivity(W/m.K)")
    plt.title("Thermal Conductivity vs Temperature for Aluminum")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
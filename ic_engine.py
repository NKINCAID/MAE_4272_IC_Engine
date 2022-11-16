"""
Simulation of a (gaseous) internal combustion engine.

"""

import cantera as ct
import numpy as np
import csv
import matplotlib.pyplot as plt

from scipy.integrate import trapz

ct.suppress_thermo_warnings()
#########################################################################
# Input Parameters
#########################################################################

# reaction mechanism, kinetics type and compositions
reaction_mechanism = "mechs/san_diego.yaml"
comp_air = "O2:1, N2:3.76"
fuel = "C3H8"

f = 1000.0 / 60.0  # engine speed [revs/s] (1000 rpm)
cr = 8.8  # compression ratio [-]
d_piston = 0.091  # piston diameter [m]
stroke = 0.086  # stroke length [m]

# Inlet temperature, pressure, and composition
T_inlet = 300.0  # K
p_inlet = ct.one_atm  # Pa
equivalence_ratio = 1.0

# outlet pressure
p_outlet = ct.one_atm  # Pa

# ambient properties
T_ambient = 300.0  # K
p_ambient = ct.one_atm  # Pa
comp_ambient = comp_air

# Inlet valve friction coefficient, open and close timings
inlet_valve_coeff = 1.0e-6
inlet_open = -18.0 / 180.0 * np.pi
inlet_close = 198.0 / 180.0 * np.pi

# Outlet valve friction coefficient, open and close timings
outlet_valve_coeff = 1.0e-6
outlet_open = 522.0 / 180 * np.pi
outlet_close = 18.0 / 180.0 * np.pi

# Fuel mass, spark open and close timings
spark_open = 353.0 / 180.0 * np.pi
spark_close = 354.0 / 180.0 * np.pi
spark_mass = 0.0  # kg

# Simulation time and parameters
sim_n_revolutions = 8
delta_T_max = 20.0
rtol = 1.0e-12
atol = 1.0e-16

#####################################################################
# Set up IC engine Parameters and Functions
#####################################################################

A_piston = 0.25 * np.pi * d_piston ** 2
V_stroke = A_piston * stroke
V_tdc = V_stroke / (cr - 1.0)


def crank_angle(t):
    """Convert time to crank angle"""
    return np.remainder(2 * np.pi * f * t, 4 * np.pi)


def piston_speed(t):
    """Approximate piston speed with sinusoidal velocity profile"""
    return -stroke / 2 * 2 * np.pi * f * np.sin(crank_angle(t))


#####################################################################
# Set up Reactor Network
#####################################################################
gas = ct.Solution(reaction_mechanism)


# define initial state and set up reactor
gas.TPX = T_inlet, p_inlet, comp_ambient
cyl = ct.IdealGasReactor(gas)
cyl.volume = V_tdc

# load reaction mechanism
gas = ct.Solution(reaction_mechanism)
# define inlet state
gas.TP = T_inlet, p_inlet
gas.set_equivalence_ratio(equivalence_ratio, fuel, comp_ambient)
inlet = ct.Reservoir(gas)

# inlet valve
inlet_valve = ct.Valve(inlet, cyl)
inlet_delta = np.mod(inlet_close - inlet_open, 4 * np.pi)
inlet_valve.valve_coeff = inlet_valve_coeff
inlet_valve.set_time_function(
    lambda t: np.mod(crank_angle(t) - inlet_open, 4 * np.pi) < inlet_delta
)

# Spark timing
spark_delta = np.mod(spark_close - spark_open, 4 * np.pi)
spark_t_open = (spark_close - spark_open) / 2.0 / np.pi / f
spark_function = lambda t: np.mod(crank_angle(t) - spark_open, 4 * np.pi) < spark_delta

# define outlet pressure (temperature and composition don't matter)
gas.TPX = T_ambient, p_outlet, comp_ambient
outlet = ct.Reservoir(gas)

# outlet valve
outlet_valve = ct.Valve(cyl, outlet)
outlet_delta = np.mod(outlet_close - outlet_open, 4 * np.pi)
outlet_valve.valve_coeff = outlet_valve_coeff
outlet_valve.set_time_function(
    lambda t: np.mod(crank_angle(t) - outlet_open, 4 * np.pi) < outlet_delta
)

# define ambient pressure (temperature and composition don't matter)
gas.TPX = T_ambient, p_ambient, comp_ambient
ambient_air = ct.Reservoir(gas)

# piston is modeled as a moving wall
piston = ct.Wall(ambient_air, cyl)
piston.area = A_piston
piston.set_velocity(piston_speed)

# create a reactor network containing the cylinder and limit advance step
sim = ct.ReactorNet([cyl])
sim.rtol, sim.atol = rtol, atol

#####################################################################
# Run Simulation
#####################################################################

# set up output data arrays
states = ct.SolutionArray(
    cyl.thermo,
    extra=("t", "ca", "V", "m", "mdot_in", "mdot_out", "dWv_dt"),
)

# simulate with a maximum resolution of 1 deg crank angle
dt = 1.0 / (360 * f)
t_stop = sim_n_revolutions / f
while sim.time < t_stop:

    if spark_function(sim.time):
        H1 = cyl.thermo.enthalpy_mass * cyl.mass

        print("Sparking now...", sim.time, crank_angle(sim.time) * 180 / np.pi)
        combustor = ct.Solution(reaction_mechanism)
        combustor.TPY = cyl.thermo.TPY
        combustor.equilibrate("UV")
        cyl.insert(combustor)
        sim.initialize()
        H2 = cyl.thermo.enthalpy_mass * cyl.mass

        sim.advance(sim.time + dt)

    else:
        # perform time integration
        sim.advance(sim.time + dt)

    # calculate results to be stored
    dWv_dt = -(cyl.thermo.P - ambient_air.thermo.P) * A_piston * piston_speed(sim.time)

    # append output data
    states.append(
        cyl.thermo.state,
        t=sim.time,
        ca=crank_angle(sim.time),
        V=cyl.volume,
        m=cyl.mass,
        mdot_in=inlet_valve.mass_flow_rate,
        mdot_out=outlet_valve.mass_flow_rate,
        dWv_dt=dWv_dt,
    )


#####################################################################
# Simple thermodynamic model
#####################################################################

V1 = V_stroke + V_tdc
P1 = ct.one_atm
T1 = 300.0
gas.TP = 300, ct.one_atm
gas.set_equivalence_ratio(equivalence_ratio, fuel, comp_air)
cp_low = gas.cp
cv_low = gas.cv
gas.equilibrate("UV")
cp_high = gas.cp
cv_high = gas.cv

cp = (cp_high + cp_low) / 2
cv = (cv_high + cv_low) / 2

gamma = cp / cv

gas.TP = 300, ct.one_atm
gas.set_equivalence_ratio(equivalence_ratio, fuel, comp_air)

V2 = V_tdc
P2 = P1 * (V1 / V2) ** gamma
T2 = T1 * (V1 / V2) ** (gamma - 1)

gas.TP = T2, P2
gas.equilibrate("UV")
V3 = V2
T3, P3 = gas.TP

V4 = V1
T4 = T1 * T3 / T2
P4 = P1 * P3 / P2

eta = 1 - 1 / (cr ** (gamma - 1))

print("Simple Thermodynamic Model:")
print("-----------------------------------------------------------")
print(
    "{:<8} | {:<10} | {:<10} | {:<10}".format("Position", "V [L]", "T [K]", "P [Bar]")
)
print("-----------------------------------------------------------")
print("{:<8} | {:<10.2f} | {:<10.2f} | {:<10.2f}".format(1, V1 * 1000, T1, P1 * 1e-5))
print("{:<8} | {:<10.2f} | {:<10.2f} | {:<10.2f}".format(2, V2 * 1000, T2, P2 * 1e-5))
print("{:<8} | {:<10.2f} | {:<10.2f} | {:<10.2f}".format(3, V3 * 1000, T3, P3 * 1e-5))
print("{:<8} | {:<10.2f} | {:<10.2f} | {:<10.2f}".format(4, V4 * 1000, T4, P4 * 1e-5))
print("-----------------------------------------------------------")
print(
    "\nEfficiency (based on gamma and compression ratio): {:.1f}%\n".format(100 * eta)
)

#######################################################################
# Plot Results in matplotlib
#######################################################################


def ca_ticks(t):
    """Helper function converts time to rounded crank angle."""
    return np.round(crank_angle(t) * 180 / np.pi, decimals=1)


t = states.t

# pressure and temperature
xticks = np.arange(0, t[-1] + 0.06, 0.06)
fig, ax = plt.subplots(nrows=2)
ax[0].plot(t, states.P / 1.0e5)
ax[0].set_ylabel("$p$ [bar]")
ax[0].set_xlabel(r"$\phi$ [deg]")
ax[0].set_xticklabels([])
ax[1].plot(t, states.T)
ax[1].set_ylabel("$T$ [K]")
ax[1].set_xlabel(r"$\phi$ [deg]")
ax[1].set_xticks(xticks)
ax[1].set_xticklabels(ca_ticks(xticks))
# plt.show()

# p-V diagram
fig, ax = plt.subplots()
ax.plot(states.V[t > 3 * t[-1] / 4] * 1000, states.P[t > 3 * t[-1] / 4] / 1.0e5)
ax.set_xlabel("$V$ [l]")
ax.set_ylabel("$p$ [bar]")
# plt.show()

# T-S diagram
fig, ax = plt.subplots()
ax.plot(
    states.m[t > 3 * t[-1] / 4] * states.s[t > 3 * t[-1] / 4],
    states.T[t > 3 * t[-1] / 4],
)
ax.set_xlabel("$S$ [J/K]")
ax.set_ylabel("$T$ [K]")
# plt.show()

# heat of reaction and expansion work
fig, ax = plt.subplots()
# ax.plot(t, 1.0e-3 * states.heat_release_rate * states.V, label=r"$\dot{Q}$")
ax.plot(t, 1.0e-3 * states.dWv_dt, label=r"$\dot{W}_v$")
# ax.set_ylim(-1e2, 1e3)
ax.legend(loc=0)
ax.set_ylabel("[kW]")
ax.set_xlabel(r"$\phi$ [deg]")
ax.set_xticks(xticks)
ax.set_xticklabels(ca_ticks(xticks))
# plt.show()

# gas composition
fig, ax = plt.subplots()
ax.plot(t, states("O2").X, label="O2")
ax.plot(t, states("CO2").X, label="CO2")
ax.plot(t, states("CO").X, label="CO")
ax.plot(t, states(fuel).X, label=fuel)
ax.legend(loc=0)
ax.set_ylabel("$X_i$ [-]")
ax.set_xlabel(r"$\phi$ [deg]")
ax.set_xticks(xticks)
ax.set_xticklabels(ca_ticks(xticks))

######################################################################
# Integral Results
######################################################################

# heat release
Q = trapz(states.heat_release_rate * states.V, t)
output_str = "{:45s}{:>4.1f} {}"
# print(
#     output_str.format(
#         "Heat release rate per cylinder (estimate):", Q / t[-1] / 1000.0, "kW"
#     )
# )

# expansion power
W = trapz(states.dWv_dt, t)
print(
    output_str.format(
        "Expansion power per cylinder (estimate):", W / t[-1] / 1000.0, "kW"
    )
)

# CO emissions
MW = states.mean_molecular_weight
CO_emission = trapz(MW * states.mdot_out * states("CO").X[:, 0], t)
CO_emission /= trapz(MW * states.mdot_out, t)

CO2_emission = trapz(MW * states.mdot_out * states("CO2").X[:, 0], t)
CO2_emission /= trapz(MW * states.mdot_out, t)


print(output_str.format("CO emission (estimate):", CO_emission * 1.0e6, "ppm"))
print(output_str.format("CO2 emission (estimate):", CO2_emission * 1.0e6, "ppm"))


# Storing data
data = []
# Corresponding data
data = np.stack(
    (
        states.t,
        states.T,
        states.P,
        states.V,
        states.s,
    ),
    axis=-1,
)
data = np.hstack((data, states.Y))

# List of data names
data_names = ["Time", "Temperature", "Pressure", "Volume", "Entropy"]
for spec_name in gas.species_names:
    data_names.append(spec_name)


filename = "data/{}_phi_{:02d}_cr_{:03d}.csv".format(
    fuel, int(equivalence_ratio * 10), int(cr * 10)
)

with open(filename, "w") as outfile:
    writer = csv.writer(outfile, delimiter=" ")
    writer.writerow(data_names)
    for n in range(len(states.t)):
        writer.writerow(data[n, :])


plt.show()

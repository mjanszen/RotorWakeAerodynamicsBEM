import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# calculate induction from the forces to get the average induction over the stream disc.

df = pd.read_csv('BEM_results.txt')
df = df.loc[199:397, :]
R_min = 10

rho = 1.225 # density of air [kg/m3]
p_static_a = np.ones(len(df))*101325 # static pressure (Point A)[N/m2]

# convert to numpy
v0 = df.v0.to_numpy()
r_b = df.r_centre.to_numpy()
r_outer = df.r_outer.to_numpy()
r_inner = df.r_inner.to_numpy()
a = df.a.to_numpy()
f_n = df.f_n.to_numpy()

omega = 1
B = 3

A_d = (r_outer**2-r_inner**2)*np.pi

# a = np.sqrt(-(f_n / (2* rho* A_d * v0**2)) + 1/4) + 1/2

# a = 1/(2*np.pi) * a


p_dynamic_a = 1/2 * rho * v0**2 # dynamic pressure (point A) [N/m2]
p_stagnation_a = p_static_a + p_dynamic_a # stagnation pressure (point A) [N/m2]
h_a = p_stagnation_a/rho + v0**2/2 # enthalpy (point A) [J = N*k = kg*m2/s2]
r_width_b = r_outer- r_inner
U_b = v0*(1-a) # velocity at disk (Point B/C)[m/s]

#%% Plot Point A
r_width_a = U_b*r_width_b / v0
r_a = np.zeros(len(r_width_a))
r_a[0] = R_min + r_width_a[0]
for i in np.arange(1,len(r_width_a)):
    r_a[i] = r_a[i-1] + r_width_a[i]
    
    
fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)

twin1 = ax.twinx()
twin2 = ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
twin2.spines.right.set_position(("axes", 1.2))


p1, = ax.plot(r_a, v0, color="black", label = "Velocity")
p2, = twin1.plot(r_a, p_stagnation_a, color="green", label = "Pressure")
p3, = twin2.plot(r_a, h_a, color="red", label = "Enthalpy")

# ax.set_xlim(10, 40)
ax.set_ylim(9.6, 10.5)
twin1.set_ylim(100000, 101500)
twin2.set_ylim(82800, 83000)

ax.set_xlabel("Radius [m]")
ax.set_ylabel("Velocity [m/s]")
twin1.set_ylabel("Pressure [N/m]")
twin2.set_ylabel("Enthalpy [J]")

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())

tkw = dict(size=4, width=1.5)
ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)

ax.grid()
# ax.legend(handles=[p1, p2, p3])
plt.savefig("figures/a_comb")
plt.show()

# plt.figure()
# plt.plot(df.v0, r_a)
# plt.ylabel("Radius [m]")
# plt.xlabel("Velocity [m/s]")
# plt.title("Velocity (Point A)")
# plt.grid()
# plt.tight_layout()
# plt.savefig("figures/v_A.png")


# plt.figure()
# plt.plot(p_stagnation_a, r_a)
# plt.ylabel("Radius [m]")
# plt.xlabel("Pressure [N/m2]")
# plt.title("Pressure (Point A)")
# plt.grid()
# plt.tight_layout()
# plt.savefig("figures/p_A.png")


#%% Plot Point B
p_stagnation_b = (h_a - U_b**2/2)*rho
h_b = h_a


fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)

twin1 = ax.twinx()
twin2 = ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
twin2.spines.right.set_position(("axes", 1.2))


p1, = ax.plot(r_b, U_b, color="black", label = "Velocity")
p2, = twin1.plot(r_b, p_stagnation_b, color="green", label = "Pressure")
p3, = twin2.plot(r_b, h_b, color="red", label = "Enthalpy")

# ax.set_xlim(0, 2)
# ax.set_ylim(0, 2)
# twin1.set_ylim(0, 4)
# twin2.set_ylim(1, 65)

ax.set_xlabel("Radius [m]")
ax.set_ylabel("Velocity [m/s]")
twin1.set_ylabel("Pressure [N/m]")
twin2.set_ylabel("Enthalpy [J]")

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())

tkw = dict(size=4, width=1.5)
ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)

ax.grid()
# ax.legend(handles=[p1, p2, p3])
plt.savefig("figures/b_comb")
plt.show()




# plt.figure()
# plt.plot(U_b, df.r_centre)
# plt.ylabel("Radius [m]")
# plt.xlabel("Velocity [m/s]")
# plt.title("Velocity (Point B)")
# plt.grid()
# plt.tight_layout()
# plt.savefig("figures/v_B.png")


# plt.figure()
# plt.plot(p_stagnation_b, df.r_centre)
# plt.ylabel("Radius [m]")
# plt.xlabel("Pressure [N/m2]")
# plt.title("Pressure (Point B)")
# plt.grid()
# plt.tight_layout()
# plt.savefig("figures/p_B.png")

#%% Plot Point C
h_c = h_b - df.f_n/rho
p_stagnation_c = (h_c - U_b**2/2)*rho
U_c = U_b

plt.figure()
plt.plot(U_c, df.r_centre)
plt.ylabel("Radius [m]")
plt.xlabel("Velocity [m/s]")
plt.title("Velocity (Point C)")
plt.grid()
plt.tight_layout()
plt.savefig("figures/v_C.png")


plt.figure()
plt.plot(p_stagnation_c, df.r_centre)
plt.ylabel("Radius [m]")
plt.xlabel("Pressure [N/m2]")
plt.title("Pressure (Point C)")
plt.grid()
plt.tight_layout()
plt.savefig("figures/p_C.png")


#%% Plot Point D
h_d = h_c
U_d = v0*(1-2*a) # velocity in far field [m/s]
p_dynamic_d = 1/2 * rho * U_d**2 # dynamic pressure (point A) [N/m2]
p_stagnation_d = p_static_a + p_dynamic_d

r_width_d = U_b*r_width_b / U_d
r_d = np.zeros(len(r_width_a))
r_d[0] = R_min + r_width_d[0]
for i in np.arange(1,len(r_width_d)):
    r_d[i] = r_d[i-1] + r_width_d[i]

plt.figure()
plt.plot(U_d, r_d)
plt.ylabel("Radius [m]")
plt.xlabel("Velocity [m/s]")
plt.title("Velocity (Point D)")
plt.grid()
plt.tight_layout()
plt.savefig("figures/v_D.png")


plt.figure()
plt.plot(p_stagnation_d, r_d)
plt.ylabel("Radius [m]")
plt.xlabel("Pressure [N/m2]")
plt.title("Pressure (Point D)")
plt.grid()
plt.tight_layout()
plt.savefig("figures/p_D.png")


#%% 
plt.figure()
plt.plot(r_a, r_a**2*np.pi, label = "Position A")
plt.plot(r_b, r_b**2*np.pi, label = "Position B")
plt.plot(r_d, r_d**2*np.pi, label = "Position D")
plt.xlabel("Radius [m]")
plt.ylabel(r"Disk Area [$m^2$]")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("figures/disk_size.png")




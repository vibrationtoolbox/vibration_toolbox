import matplotlib.pyplot as plt
import vibration_toolbox as vtb
omega_n, x, U = vtb.euler_beam_modes(n=1)
plt.plot(x,U)
plt.xlabel('x (m)')
plt.ylabel('Displacement (m)')
plt.title('Mode 1')

import matplotlib.pyplot as plt
import vibration_toolbox as vtb
omega_n, x, U = vtb.euler_beam_modes(n=1)
plt.figure()
# <matplotlib.figure...>
plt.plot(x,U)
# [<matplotlib.lines.Line2D object at ...>]
plt.xlabel('x (m)')
# Text(0.5,0,'x (m)')
plt.ylabel('Displacement (m)')
# Text(0,0.5,'Displacement (m)')
plt.title('Mode 1')
# Text(0.5,1,'Mode 1')
plt.grid('on')

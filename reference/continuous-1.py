import matplotlib.pyplot as plt
import vibration_toolbox as vtb
omega_n, x, U = vtb.euler_beam_modes(n=1)
plt.figure()
# <matplotlib.figure...>
plt.plot(x,U)
# [<matplotlib.lines.Line2D object at ...>]
plt.xlabel('x (m)')
# <matplotlib.text.Text object at ...>
plt.ylabel('Displacement (m)')
# <matplotlib.text.Text object at ...>
plt.title('Mode 1')
# <matplotlib.text.Text object at ...>
plt.grid('on')

import matplotlib.pyplot as plt
import vibration_toolbox as vtb
vtb.free_response()[1][:5] # get the first five values of x
# array([[ 1.  ],
# [ 1.  ],
# [ 0.99],
# [ 0.99],
# [ 0.98]])

t, x, *_ = vtb.free_response() # *_ ignores all other returns
plt.plot(t,x)
# [<matplotlib.lines.Line2D object at ...>]
plt.xlabel('Time (sec)')
# <matplotlib.text.Text object at ...>
plt.ylabel('Displacement (m)')
# <matplotlib.text.Text object at ...>
plt.title('Displacement versus time')
# <matplotlib.text.Text object at ...>
plt.grid('on')

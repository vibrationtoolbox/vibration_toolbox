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
# Text(0.5,0,'Time (sec)')
plt.ylabel('Displacement (m)')
# Text(0,0.5,'Displacement (m)')
plt.title('Displacement versus time')
# Text(0.5,1,'Displacement versus time')
plt.grid('on')

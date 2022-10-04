from DSP import DSP
import numpy as np
x = np.array(range(0,1000))
y = list(x)
for i in list(x):
    y[i] = 0.1*np.sin((2*np.pi*400/1000)*i)+1*np.sin((2*np.pi*10/1000)*i) + 0.5

filt = (0.5/0.537)*np.array([
  -0.02062630863958946,
  -0.010688729607911805,
  0.09988368103096376,
  0.2805270356057702,
  0.3730973942070046,
  0.2805270356057702,
  0.09988368103096376,
  -0.010688729607911805,
  -0.02062630863958946
])
DSP = DSP(data = y, fs = 1000, filt = 'Yes', num = filt)
DSP.calculate()
DSP.plot()
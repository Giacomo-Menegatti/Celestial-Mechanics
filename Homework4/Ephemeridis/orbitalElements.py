import numpy as np
class body:
  pass

mercury = body()
mercury.name = 'Mercury'
mercury.L = np.array([252.250906, 149474.0722491, 0.00030350, 0.000000018])
mercury.a = np.array([0.387098310,     0.0,          0.0,          0.0])
mercury.e = np.array([0.20563175,      +0.000020407,  -0.0000000283, -0.00000000018])
mercury.i = np.array([7.004986,      +0.0018215,    -0.00001810,   +0.000000056])
mercury.Omega = np.array([48.330893,      +1.1861883,    +0.00017542,  +0.000000215])
mercury.pi = np.array([77.456119   ,      +1.5564776,    +0.00029544,   +0.000000009])

venus = body()
venus.name = 'Venus'
venus.L = np.array([181.979801,     +58519.2130302  ,  +0.00031014 ,  +0.000000015 ])
venus.a = np.array([  0.723329820,      0.0         ,  0.0         ,  0.0          ])
venus.e = np.array([  0.00677192 ,      -0.000047765,  +0.0000000981, +0.00000000046])
venus.i = np.array([  3.394662   ,      +0.0010037  ,  -0.00000088 ,  -0.000000007  ])
venus.Omega = np.array([ 76.679920   ,      +0.9011206  ,  +0.00040618 , -0.000000093   ])
venus.pi = np.array([ 131.563703  ,       +1.4022288 ,   -0.00107618 ,  -0.000005678  ])

earth=body()
earth.name = 'Earth'
earth.L= np.array([100.466457   ,  +36000.7698278   , +0.00030322   ,+0.000000020     ])
earth.a= np.array([  1.000001018,       0.0         ,  0.0          , 0.0             ])
earth.e= np.array([  0.01670863 ,      -0.000042037 , -0.0000001267 ,+0.00000000014   ])
earth.i= np.array([  0.0        ,       0.0         ,  0.0          , 0.0             ])
earth.Omega= np.array([  0.0        ,       0.0         ,  0.0          , 0.0             ])
earth.pi  = np.array([102.937348   ,      +1.7195366   , +0.00045688   ,-0.000000018     ])

mars=body()
mars.name = 'Mars'
mars.L = np.array([355.433000   ,  +19141.6964471    ,+0.00031052   ,+0.000000016     ])
mars.a = np.array([ 1.523679432,       0.0          , 0.0          , 0.0             ])
mars.e = np.array([ 0.09340065 ,      +0.000090484  ,-0.0000000806 ,-0.00000000025   ])
mars.i = np.array([ 1.849726   ,      -0.0006011    ,+0.00001276   ,-0.000000007     ])
mars.Omega = np.array([49.558093   ,      +0.07720959   ,+0.00001557   ,+0.00002267      ])
mars.pi = np.array([336.060234   ,      +1.8410449    ,+0.00013477   ,+0.000000536     ])

sun = body()
sun.name = 'sun'
sun.L = np.array([0,0,0,0])
sun.a = np.array([0,0,0,0])
sun.e = np.array([1,0,0,0])
sun.i = np.array([0,0,0,0])
sun.Omega = np.array([0,0,0,0])
sun.pi = np.array([0,0,0,0])
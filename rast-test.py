from scipy import *
from numpy import *
import rast

lines = rast.makeLinesP2D()
lines.set_maxresults(5)
lines.set_error(2.0,10.0)
lines.set_tolerance(0.1,0.0001)
lines.set_minweight(4)
lines.set_maxoffset(250)
lines.set_lsq(1)

for i in range(20):
    x = random.uniform(0.0,100.0)
    y = random.uniform(0.0,100.0)
    lines.add_ipoint(x,y,0.0,1.0)
for i in range(20):
    x = random.uniform(0.0,100.0)
    y = 20.0 + x * 0.3
    lines.add_ipoint(x,y,0.0,1.0)

lines.compute()

for i in range(lines.nresults()):
    print lines.angle(i),lines.offset(i),lines.weight(i),lines.nmatches(i)

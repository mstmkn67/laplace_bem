from numpy import *
from UDFManager import *
import matplotlib.pyplot as plt

number=range(1,11)
number_face=[]
calc_time0=[]
calc_time1=[]
for i in number:
	udf=UDFManager("./line%02d_o.udf"%(i))
	f=udf.size("input.face[]")
	t=udf.get("output.time")
	number_face.append(f)
	calc_time0.append(t)
for i in number:
	udf=UDFManager("./line%02d_fo.udf"%(i))
	f=udf.size("input.face[]")
	t=udf.get("output.time")
	calc_time1.append(t)
fig = plt.figure(figsize=(12.0, 9.0))
ax=fig.add_subplot()
ax.plot(number_face,calc_time0,marker='o',markersize=5,linestyle='-',color='blue')
ax.plot(number_face,calc_time1,marker='o',markersize=5,linestyle='-',color='red')
ax.set_xlabel("number of faces")
ax.set_ylabel("calculation time")
plt.tight_layout()
plt.savefig("calculation_time.png")
plt.show()


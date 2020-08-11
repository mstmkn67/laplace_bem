from numpy import *
from UDFManager import *
import matplotlib.pyplot as plt
def outer(a,b):
	return array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])

number=range(1,11)
evaporation_rate=[]
for i in number:
	udf=UDFManager("./line%02d_fo.udf"%(i))
	v=array(udf.get("input.vertex[].position"))
	f=udf.get("input.face[].vertex[]")
	q=udf.get("output.face[].q")
	D=26.10
	a=0.0
	for j in range(len(f)):
		r0,r1,r2=v[f[j][0]],v[f[j][1]],v[f[j][2]]
		s=outer(r1-r0,r2-r0)
		s=0.5*sqrt(dot(s,s))
		a+=s*q[j]
	evaporation_rate.append(-D*a/2.0/len(udf.get("input.region_condition[]")))
fig = plt.figure(figsize=(12.0, 9.0))
ax=fig.add_subplot()
ax.plot(number,evaporation_rate,marker='o',markersize=5,linestyle='-',color='red')
ax.set_xlabel("droplet number")
ax.set_ylabel("evaporation rate per a droplet [ug]")
plt.tight_layout()
plt.savefig("evaporation_rate.png")
plt.show()


from numpy import *
def outer(a,b):
	return array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
q=$output.face[].q
face=$input.face[].vertex[]
for i in range(len(q)):
	r0=array(get(getLocation("Vertex",face[i][0])+".position"))
	r1=array(get(getLocation("Vertex",face[i][1])+".position"))
	r2=array(get(getLocation("Vertex",face[i][2])+".position"))
	c=(r0+r1+r2)/3.0
	n=outer(r1-r0,r2-r0)
	n/=sqrt(dot(n,n))
	line(c.tolist(),(c-5*q[i]*n).tolist(),1)
	polyline([r0.tolist(),r1.tolist(),r2.tolist(),r0.tolist()],[1,0,0,0.5])

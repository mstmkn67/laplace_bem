$input.Gauss_Legendre_integral_points_number="7"
$input.bicgstab.tolerance=1e-5
import MMesh5
from numpy import *
cv,H=2.32e-2,0.4 
$input.external_field.phi=H*cv
$input.external_field.gradient_phi=[0.0,0.0,0.0]

N=1
DX=3.0#mm
pos,angle,radius=[],[],[]
for i in range(N):
	pos.append([i*DX,0.0,0.0])
	angle.append(41.0*pi/180.)
	radius.append(1.0)#mm
$input.system.center=[DX*(N-1)*0.5,0.0,0.0]
$input.system.size=DX*N
$input.vertex[]=[]
$input.face[]=[]
$input.region_condition[]=[]
$input.evaluation_point[]=[]


#vertex,face=MMesh5.sphere(1,3,"octahedron")
vertex1,face1=MMesh5.disk(1,4,"hexagon")
vertex2,face2=MMesh5.disk(1,4,"hexagon")
for i in range(len(vertex1)):
	vertex1[i][0],vertex1[i][1]=vertex1[i][1],vertex1[i][0]
	vertex1[i][2]=1.0
	vertex2[i][2]=-1.0
n1=len(vertex1)
for i in range(len(face2)):
	face2[i][0]+=n1
	face2[i][1]+=n1
	face2[i][2]+=n1
vertex=array(vertex1.tolist()+vertex2.tolist())
face=face1+face2

nv=len(vertex)
nf=len(face)
for n in range(N):
	h0=radius[n]*tan(0.5*angle[n])
	for i in range(nv):
		$input.vertex[n*nv+i].id=n*nv+i
		r2=vertex[i][0]*vertex[i][0]+vertex[i][1]*vertex[i][1]
		s2=sin(angle[n])*sin(angle[n])
		z=sqrt(radius[n]*radius[n]/s2-r2)-radius[n]/tan(angle[n])
		if vertex[i][2]!=0.0:
			vertex[i][2]*=z/abs(vertex[i][2])
		$input.vertex[n*nv+i].position=(radius[n]*vertex[i]+pos[n]).tolist()
	$input.region_condition[n].name="%d"%(n)
	$input.region_condition[n].type="phi"
	$input.region_condition[n].phi=cv
	for i in range(nf):
		$input.face[nf*n+i].id=nf*n+i
		$input.face[nf*n+i].vertex[]=[nv*n+face[i][0],nv*n+face[i][1],nv*n+face[i][2]]
		$input.region_condition[n].face[i]=n*nf+i



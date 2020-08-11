
import random
import MMesh5
from numpy import *
N=10
rr=1.0
sep=2.1
ss=10.0
pos=[]
for i in range(N):
	while(1):
		flag=True
		p=array([ss*(random.random()-0.5),ss*(random.random()-0.5),ss*(random.random()-0.5)])
		for j in pos:
			if dot(j-p,j-p)<sep*sep:
				flag=False
				break
		if flag==True:
			break
	pos.append(p)
#for i in pos:
#	sphere(i.tolist(),[1,0,0,1,rr])
#for i in range(len(pos)):
#	for j in range(i+1,len(pos)):
#		print dot(pos[i]-pos[j],pos[i]-pos[j])

$input.vertex[]=[]
$input.face[]=[]
$input.region_condition[]=[]
$input.evaluation_point[]=[]

vertex,face=MMesh5.sphere(1,3,"icosahedron")
nv=len(vertex)
nf=len(face)
for n in range(N):
	for i in range(nv):
		$input.vertex[n*nv+i].id=n*nv+i
		$input.vertex[n*nv+i].position=(rr*vertex[i]+pos[n]).tolist()
	$input.region_condition[n].name="%d"%(n)
	$input.region_condition[n].type="q"
	$input.region_condition[n].q=0
	for i in range(nf):
		$input.face[nf*n+i].id=nf*n+i
		$input.face[nf*n+i].vertex[]=[nv*n+face[i][0],nv*n+face[i][1],nv*n+face[i][2]]
		$input.region_condition[n].face[i]=n*nf+i

ss=7.0
vertex,face=MMesh5.rect3d(ss,ss,ss,4)
n=len(vertex)
for i in range(n):
	$input.vertex[N*nv+i].id=N*nv+i
	$input.vertex[N*nv+i].position=vertex[i].tolist()
m=len(face)
for i in range(m):
	$input.face[N*nf+i].id=N*nf+i
	$input.face[N*nf+i].vertex[]=[N*nv+face[i][0],N*nv+face[i][2],N*nv+face[i][1]]
region_xym,region_xyp,region_yzm,region_yzp,region_zxm,region_zxp=[],[],[],[],[],[]
for i in range(m):
	r0,r1,r2=vertex[face[i][0]],vertex[face[i][1]],vertex[face[i][2]]
	c=(r0+r1+r2)/3.0
	if c[2]==-ss: #xy-
		region_xym.append(nf*N+i)
	if c[2]==ss: #xy+
		region_xyp.append(nf*N+i)
	if c[0]==-ss: #yz-
		region_yzm.append(nf*N+i)
	if c[0]==ss: #yz+
		region_yzp.append(nf*N+i)
	if c[1]==-ss: #zx-
		region_zxm.append(nf*N+i)
	if c[1]==ss: #zx+
		region_zxp.append(nf*N+i)
$input.region_condition[N+0].name="zm"
$input.region_condition[N+0].face[]=region_xym
$input.region_condition[N+0].type="q"
$input.region_condition[N+0].q=0.0
#
$input.region_condition[N+1].name="zp"
$input.region_condition[N+1].face[]=region_xyp
$input.region_condition[N+1].type="q"
$input.region_condition[N+1].q=0.0
#
$input.region_condition[N+2].name="xm"
$input.region_condition[N+2].face[]=region_yzm
$input.region_condition[N+2].type="phi"
$input.region_condition[N+2].phi=0.0
#
$input.region_condition[N+3].name="xp"
$input.region_condition[N+3].face[]=region_yzp
$input.region_condition[N+3].type="phi"
$input.region_condition[N+3].phi=100.0
#
$input.region_condition[N+4].name="ym"
$input.region_condition[N+4].face[]=region_zxm
$input.region_condition[N+4].type="q"
$input.region_condition[N+4].q=0.0
#
$input.region_condition[N+5].name="yp"
$input.region_condition[N+5].face[]=region_zxp
$input.region_condition[N+5].type="q"
$input.region_condition[N+5].q=0.0


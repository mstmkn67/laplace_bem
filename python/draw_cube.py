from numpy import *
def draw_cube_polygon(center,size,color,udf):
	r=array([[ 0.5, 0.5, 0.5],[-0.5, 0.5, 0.5],[-0.5,-0.5, 0.5],[ 0.5,-0.5, 0.5],
             [ 0.5, 0.5,-0.5],[-0.5, 0.5,-0.5],[-0.5,-0.5,-0.5],[ 0.5,-0.5,-0.5]])
	for i in range(8):
		r[i]=size*r[i]+center
	r=r.tolist()
	udf.polygon([r[0],r[1],r[2],r[3]],color)
	udf.polygon([r[4],r[7],r[6],r[5]],color)
	udf.polygon([r[0],r[3],r[7],r[4]],color)
	udf.polygon([r[3],r[2],r[6],r[7]],color)
	udf.polygon([r[1],r[2],r[6],r[5]],color)
	udf.polygon([r[1],r[0],r[4],r[5]],color)

def draw_cube_polyline(center,size,color,udf):
	r=array([[ 0.5, 0.5, 0.5],[-0.5, 0.5, 0.5],[-0.5,-0.5, 0.5],[ 0.5,-0.5, 0.5],
             [ 0.5, 0.5,-0.5],[-0.5, 0.5,-0.5],[-0.5,-0.5,-0.5],[ 0.5,-0.5,-0.5]])
	for i in range(8):
		r[i]=size*r[i]+center
	r=r.tolist()
	udf.polyline([r[0],r[1],r[2],r[3],r[0]],color)
	udf.polyline([r[4],r[7],r[6],r[5],r[4]],color)
	udf.polyline([r[0],r[3],r[7],r[4],r[0]],color)
	udf.polyline([r[3],r[2],r[6],r[7],r[3]],color)
	udf.polyline([r[1],r[2],r[6],r[5],r[1]],color)
	udf.polyline([r[1],r[0],r[4],r[5],r[1]],color)

#draw_cubic_polygon(array([0,0,0]),1,1)
#draw_cubic_polyline(array([0,1,2]),1,2)

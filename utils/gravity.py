from Dynearthsol import Dynearthsol
import matplotlib.pyplot as plt 
import matplotlib.tri as tri
import numpy as np 
import triangle as tr
import sys 

G=6.67e-11
modelname="RC17"
modelinfo=modelname+".info"
nframe=9999
height=3000 
left=30000
right=110000

def area(Ecoord):
	ab=Ecoord[:,1,:]-Ecoord[:,0,:]
	ac=Ecoord[:,2,:]-Ecoord[:,0,:]
	_vec =np.hstack([ab,ac]).reshape(-1,2,2)
	area=np.abs(np.linalg.det(_vec))/2
	return area


def gravity(pij,Ecenter,area,rho):
	dv=pij-Ecenter
	dr=np.sum(dv**2,axis=1)
	g0=2*G*np.dot(dv[:,1]*area/dr,rho)
	return g0

def freeair(vec_pij,Ecoord,rho):
	nsmaple=vec_pij.shape[0]
	g=np.empty(shape=nsmaple)
	A0=area(Ecoord)
	Ecenter=np.average(Ecoord,axis=1)
	for p0 in range(nsmaple):
		g[p0]=gravity(vec_pij[p0],Ecenter,A0,rho)
	return g

def bouguer(vx,vy,coord1,conn,rho):
	Ecoord=coord1[conn]
	nsmaple=len(vx)
	g=np.empty(shape=nsmaple)
	Ecenter=np.average(Ecoord,axis=1)
	mask=Ecenter[:,1]<0
	Ecenter=Ecenter[mask]
	Ecoord=Ecoord[mask]
	rho=rho[mask]
	A0=area(Ecoord)
	for p0 in range(nsmaple):
		g[p0]=gravity([vx[p0],vy[p0]],Ecenter,A0,rho)
	# g=g+0.386*vec_pij[:,1]
	return g

# def talwani(coord,conn,pij,rho):
# 	Ecoord= coord[conn]
# 	phi1=Ecoord[:,1,:]-Ecoord[:,0,:]
# 	phi2=Ecoord[:,2,:]-Ecoord[:,1,:]
# 	phi3=Dcoord[:,0,:]-Ecoord[:,2,:]
# 	dij=Ecoord[:,0,:]-pij
# 	theta=np.tan(dij[:,1]/dij[:,0])
# 	dg=2*G*rho
# 	return 


def main(): 
	mopen=Dynearthsol(modelname)
	rho=mopen.read_field(nframe,'density')
	conn=mopen.read_field(nframe,'connectivity')
	coord=mopen.read_field(nframe,'coordinate')
	bcflag=mopen.read_field(nframe,'bcflag').astype(np.int16)
	node32=coord[(bcflag&32).astype(bool)]



	###### create new vertice #######
	## extend two walls 
	xbtm=coord[bcflag==17]
	xbtm2=[[xbtm[0,0]-10000,xbtm[0,1]]]
	xtop=coord[bcflag==33]
	xtop2=[[xtop[0,0]-10000,xtop[0,1]]]
	x0=np.vstack([xbtm,xbtm2,xtop, xtop2])
	tridict=dict(vertices=x0)
	B=tr.triangulate(tridict,'qa1000000')
	n0=coord.shape[0]
	coord=np.vstack([coord,B['vertices']])
	conn=np.vstack([conn,B['triangles']+n0])	
	rho=np.hstack([rho,np.ones(shape=B['triangles'].shape[0])*1750])
	bcflag=np.hstack([bcflag,np.zeros(shape=B['vertices'].shape[0])]).astype(np.int16)

	xbtm=coord[bcflag==66]
	xbtm2=[[xbtm[0,0]+10000,xbtm[0,1]]]
	xtop=coord[bcflag==34]
	xtop2=[[xtop[0,0]+10000,xtop[0,1]]]
	x0=np.vstack([xbtm,xbtm2,xtop, xtop2])
	tridict=dict(vertices=x0)
	B=tr.triangulate(tridict,'qa1000000')
	n0=coord.shape[0]
	coord=np.vstack([coord,B['vertices']])
	conn=np.vstack([conn,B['triangles']+n0])
	rho=np.hstack([rho,np.ones(shape=B['triangles'].shape[0])*1750])

	################################

	
	measureX=np.arange(left,right,500)
	measureY=np.ones(shape=measureX.shape)*height



	result=bouguer(measureX,measureY,coord,conn,rho)*1e5



	##### figure #####
	fig,ax=plt.subplots(figsize=(10,2.5))
	ax.set_title(modelname)
	ax.plot(measureX,result)
	ax.set_ylabel('Bouguer gravity anomaly (mgal)(blue)')
	ax2=ax.twinx()
	ax2.plot(node32[:,0],node32[:,1],'r')
	ax2.set_ylabel('topography (m)(red)')
	# ax2.set_aspect('equal', adjustable="datalim")
	plt.show()

	
if __name__== "__main__":
	main()
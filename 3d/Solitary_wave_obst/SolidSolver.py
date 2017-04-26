#!/usr/bin/env python3
#Filename: SolidSolver.py
from dolfin import *
import mshr
import numpy as np
beta, gamma, eta = 0.25, 0.5, 0.2
class ElasticBody():
    def __init__(self, E=200e9, nu=0.3, g=np.array([0., 0., -9.8]), density= 7800.):
        self.E=E
        self.nu=nu
        self.g=np.array([0., 0., -9.8 ])
        self.density=density
        self.mu=E/(2.*(1.+nu))
        self.lmbda=E*nu/((1. + nu)*(1. - 2.*nu))
    def sigma(self, u):#define the stress
        return 2.0*self.mu*sym(grad(u)) + self.lmbda*tr(grad(u))*Identity(u.geometric_dimension())
 
class ElastoDynBody(ElasticBody):
    def __init__(self, E=200e9, nu=0.3, g=np.array([0., 0., -9.8]), density= 7800.,eta=0.2):
        super(ElastoDynBody, self).__init__(E, nu, g, density)
        self.eta=eta
        #self.E=E
        #self.nu=nu
        #self.g=np.array([0., 0., -9.8 ])
        #self.density=density
		#self.mu=E/(2.*(1.+nu))
        #self.lmbda=E*nu/((1. + nu)*(1. - 2.*nu))
    def sigma(self, u, v):#define the stress
        return 2.0*self.mu*sym(grad(u)) + (self.lmbda*tr(grad(u)) + self.eta*tr(grad(v)))*Identity(u.geometric_dimension())
         
class ElastoDynCylinder(ElastoDynBody):
    def __init__(self,E=200e9,nu=0.33,g=np.array([0.,0.,-9.8]),density=7800,D=1.,H=1.,loc=np.array([10,10])):
		super(ElastoDynCylinder, self).__init__(E, nu, g, density)
		self.diameter=D
		self.H=H
		self.loc=loc


    def GenerateMesh(self,he):
        mesh2D = mshr.Circle(Point(self.loc[0], self.loc[1], 0), self.D/2.)
        mesh3D = mshr.Extrude2D(mesh2D, self.H)
        self.mesh= mshr.generate_mesh(mesh3D, int(np.pi*self.D/he)+1)
	def set_dt(self, dt):
		self.dt=dt

	def set_functional(self, FuncType='Lagrange', Order=1):
		self.V=VectorFunctionSpace(self.mesh, FuncType, Order)
        self.u1, self.w = TrialFunction(V), TestFunction(V)
        self.u0, self.v0, self.a0 = Function(V), Function(V),Function(V)
        self.v1 = (gamma/(beta*self.dt))*(self.u1 - self.u0) - (gamma/beta - 1.0)*self.v0 - self.dt*(gamma/(2.0*beta) - 1.0)*self.a0
        self.a1 = (1.0/(beta*self.dt**2))*(self.u1 - self.u0 - self.dt*self.v0) - (1.0/(2.0*beta) - 1.0)*self.a0
    def set_force(force_x=None,force_y=None,force_z=None,force_coords=None):
        assert force_x.shape[0]==force_y.shape[0]==force_z.shape[0]==force_coords.shape[0]
        assert force_x.shape[1]==force_y.shape[1]==force_z.shape[1]==force_coords.shape[1]
        self.force_x=force_x
        self.force_y=force_y
        self.force_z=force_z
        self.force_coords=force_coords
        
    def set_variationalform():
        h=self._traction() 
        self.F=(self.density*dot(self.a1, self.w)+inner(sigma(self.u1, self.v1),grad(self.w)))*dx-dot(h, self.w)*ds
        self.a,self.l=lhs(self.F), rhs(self.F)

	def set_DBC(BC={}):
        self.DBC=[]
		for boundary, value in BC.iteritems():
			boundaries=lambda x, on_boundary: on_boundary and boundary
            bc=DirichletBC(self.V, Constant(value[0],value[1],value[2]),boundaries)
			self.DBC.append(bc)
    def solve():
        self.uk=Function(self.V)
        solve(self.a==self.l,self.uk,self.DBC)





    
    def _traction(Expression):
        def __init__(self):
            pass
        def eval(self,value,x):
            from scipy import spatial
            nObstEleBoundary=force_coords.shape[0]
            nEleBoundaryQPoint=force_coords.shape[1]
            Dim=force_coords.shape[2]
            
            tree=spatial.KDTree(np.reshape(self.force_coords,(nObstEleBoundary*nEleBoundaryQPoint,Dim)))
            value=self.force_x[tree.query(x)[1]] #nearest neighbor
            #value[1]=
            #value[2]=

        def value_shape(self):
            return (3,)



class PlasticCylinder:
	pass

class HyperelastoDynCylinder(ElastoDynCylinder):
	pass    

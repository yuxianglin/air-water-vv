#!/usr/bin/env python3
#Filename: SolidSolver.py
from dolfin import *
import mshr
import numpy as np
from proteus.Profiling import logEvent
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
class ElasticBody(object):
    def __init__(self, E=200e9, nu=0.3, g=np.array([0., 0., -9.8]), density= 7800.):
        self.E=E
        self.nu=nu
        self.g=g
        self.density=density
        self.mu=E/(2.*(1.+nu))
        self.lmbda=E*nu/((1. + nu)*(1. - 2.*nu))

    def _sigma(self, u):#define the stress
        return 2.0*self.mu*sym(grad(u)) + self.lmbda*tr(grad(u))*Identity(u.geometric_dimension())

    def set_force(self,force_x=None,force_y=None,force_z=None,force_coords=None):
        assert force_x.shape[0]==force_y.shape[0]==force_z.shape[0]==force_coords.shape[0]
        assert force_x.shape[1]==force_y.shape[1]==force_z.shape[1]==force_coords.shape[1]
        self.force_x=force_x
        self.force_y=force_y
        self.force_z=force_z
        self.force_coords=force_coords

    def set_DBC(self,BC={}):
        self.DBC=[]
        for boundary, value in BC.iteritems():
            #boundaries=boundary
            bc=DirichletBC(self.V, Constant((value[0],value[1],value[2])), boundary)
            self.DBC.append(bc)

    def set_variationalform(self):
        h=_traction(self.force_coords,self.force_x,self.force_y,self.force_z)
        f=Constant((self.g[0],self.g[1],self.g[2]))
        self.F=inner(self._sigma(self.u1),grad(self.w))*dx-dot(h,self.w)*ds-self.density*dot(f,self.w)*dx
        self.a,self.l=lhs(self.F), rhs(self.F)
        #self.a=inner(self._sigma(self.u1),grad(self.w))*dx
        #self.l=dot(h,self.w)*ds+self.density*dot(f,self.w)*dx

    def Gen_Tfunc(self):
        self.u1, self.w = TrialFunction(self.V), TestFunction(self.V)
        self.u = Function(self.V)

    def solve(self):
        solve(self.a==self.l,self.u,self.DBC)
        file=File('u.pvd')
        file << self.u
        logEvent('displacement write to u.pvd')

class ElasicCylinder(ElasticBody):
    def __init__(self,E=200e9, nu=0.3, g=np.array([0., 0., -9.8]), density= 7800.,D=1.,H=1.,loc=np.array([10,10])):
        super(ElasicCylinder,self).__init__(E, nu, g, density)
        self.diameter=D
        self.H=H
        self.loc=loc

    def Gen_Mesh(self,he=0.05,FuncType='Lagrange', Order=1,output=0,mesh=None):
        if mesh:
            self.mesh=mesh
        else:
            mesh2D = mshr.Circle(Point(self.loc[0], self.loc[1], 0), self.diameter/2.)
            mesh3D = mshr.Extrude2D(mesh2D, self.H)
            self.mesh= mshr.generate_mesh(mesh3D, int(np.pi*self.diameter/he)+1)
        self.V=VectorFunctionSpace(self.mesh, FuncType, Order)
        if output:
            file=File('mesh.pvd')
            file << self.mesh
            logEvent("mesh.pvd generated")




class ElastoDynBody(ElasticBody):
    def __init__(self, E=200e3, nu=0.3, g=np.array([0., 0., -9.8]), density= 7800.,beta =0.25, gamma = 0.5):
        super(ElastoDynBody, self).__init__(E, nu, g, density)
        self.eta=self.lmbda/20.
        self.beta=beta
        self.gamma=gamma

    def _sigma(self, u, v):#define the stress
        """define the stress considering the damping"""
        return 2.0*self.mu*sym(grad(u)) + (self.lmbda*tr(grad(u))+ self.eta*tr(grad(v)))*Identity(u.geometric_dimension())

    def set_function(self):
        pass

    def set_dt(self, dt):
        self.dt=dt

    def set_variationalform(self):
        self.v1 = (self.gamma/(self.beta*self.dt))*(self.u1 - self.u0)- (self.gamma/self.beta - 1.0)*self.v0- self.dt*(self.gamma/(2.0*self.beta) - 1.0)*self.a0
        self.a1 = (1.0/(self.beta*self.dt**2))*(self.u1 - self.u0 - self.dt*self.v0)- (1.0/(2.0*self.beta) - 1.0)*self.a0

        h=_traction(self.force_coords,self.force_x,self.force_y,self.force_z)
        f=Constant((self.g[0],self.g[1],self.g[2]))
        self.F=(self.density*dot(self.a1, self.w)+inner(self._sigma(self.u1, self.v1),grad(self.w)))*dx- dot(h, self.w)*ds -self.density*dot(f,self.w)*dx
        self.a,self.l=lhs(self.F), rhs(self.F)

    def _update(self):
        u_vec, u0_vec = self.u.vector(),self.u0.vector()
        v0_vec, a0_vec = self.v0.vector(),self.a0.vector()
        a_vec = (1.0/(2.0*self.beta))*((u_vec - u0_vec -v0_vec*self.dt)/(0.5*self.dt**2) - (1.0-2.0*self.beta)*a0_vec)
        v_vec = self.dt*((1.0-self.gamma)*a0_vec + self.gamma*a_vec) + v0_vec
        #updata tn-->tn+1
        self.v0.vector()[:], self.a0.vector()[:] = v_vec, a_vec
        self.u0.vector()[:]=u_vec

    def solve(self):
        solve(self.a==self.l,self.u,self.DBC)
        self._update()
        file=File('u1.pvd')
        file << self.u
        logEvent('displacement write to u.pvd')

    def Gen_Tfunc(self):
        self.u1, self.w = TrialFunction(self.V), TestFunction(self.V)
        self.u = Function(self.V)
        self.u0, self.v0, self.a0 = Function(self.V), Function(self.V),Function(self.V)


    def set_Init(self,u=Constant((0.0,0.0,0.0)),v=0):
        self.u0.vector()[:]=u.vector()[:]# interpolate u
        self.v0.vector().array()[:]=v# hard co


class ElastoDynCylinder(ElastoDynBody):
    def __init__(self,E=200e9,nu=0.33,g=np.array([0.,0.,-9.8]),density=7800,D=1.,thick=0,H=1.,loc=np.array([10,10])):
        super(ElastoDynCylinder, self).__init__(E, nu, g, density)
        self.diameter=D
        self.H=H
        self.thick=thick
        self.loc=loc
        
    def Gen_Mesh(self,he=0.05,FuncType='Lagrange', Order=1,output=0,mesh=None):
        if mesh:
            self.mesh=mesh
        else:
            mesh2D= mshr.Circle(Point(self.loc[0], self.loc[1], 0), self.diameter/2.)
            if self.thick:
                mesh2D=mesh2D-mshr.Circle(Point(self.loc[0], self.loc[1], 0), self.diameter/2.-self.thick)
            mesh3D = mshr.Extrude2D(mesh2D, self.H)
            self.mesh= mshr.generate_mesh(mesh3D, int(np.pi*self.diameter/he)+1)
        self.V=VectorFunctionSpace(self.mesh, FuncType, Order)
        if output:
            file=File('mesh.pvd')
            file << self.mesh
            logEvent("mesh.pvd generated")

class _traction(Expression):
    def __init__(self,force_coords,force_x,force_y,force_z):
        #super(_traction,self).__init__(self)
        from scipy import spatial
        self.force_coords=force_coords
        nObstEleBoundary=force_coords.shape[0]
        nEleBoundaryQPoint=force_coords.shape[1]
        Dim=force_coords.shape[2]
        self.force_x=np.reshape(force_x,nObstEleBoundary*nEleBoundaryQPoint)
        self.force_y=np.reshape(force_y,nObstEleBoundary*nEleBoundaryQPoint)
        self.force_z=np.reshape(force_z,nObstEleBoundary*nEleBoundaryQPoint)
        self.tree=spatial.KDTree(np.reshape(self.force_coords,(nObstEleBoundary*nEleBoundaryQPoint,Dim)))

    def eval(self,value,x):
        value[:]=0

        ind=self.tree.query(x)[1]
        d=self.tree.query(x)[0]
        value[0]=self.force_x[ind] #nearest neighbor
        value[1]=self.force_y[ind] #nearest neighbor
        value[2]=self.force_z[ind] #nearest neighbor

    def value_shape(self):
        return (3,)



class PlasticCylinder:
	pass

class HyperelastoDynCylinder(ElastoDynCylinder):
	pass    

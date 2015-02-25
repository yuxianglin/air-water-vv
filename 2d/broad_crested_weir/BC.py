class boundaryConditions:

# Checking for entering the correct BC type
    def BCTypeCheck(self,BCType):
        boundaryList = ["pDirichlet","uDirichlet","vDirichlet", "wDirichlet","pAdvective","uAdvective","vAdvective","wAdvective","uDiffusive","vDiffusive","wDiffusive","vofDirichlet","vofAdvective"]
        if BCType not in boundaryList:
            print("Boundary condition type not valid")
            exit(1)


# Basic functions
    def empty(self):
        return None    
    def constantValue(self,value):
        return lambda x,t: value
    def linear(self,a0,a1,i):
        return lambda x,t: a0 + a1*x[i]
    def phase(self,k,ww,cos=True):
        if cos:
            return lambda x,t: cos(k*x-ww*t)
        else:
            return lambda x,t: sin(k*x-ww*t)
    def hyperPhase(self,k,ww,cosh=True):
        if cosh:
            return lambda x,t: cosh(k*x-ww*t)
        else:
            return lambda x,t: sinh(k*x-ww*t)


#Basic boundary conditions
#NoSlip: U=0
    def freeSlip(self,BCType):
        """Imposes a free slip condition
        Takes the following arguments
        BCType: Type of boundary condition 
        Returns Advective conditions as zero, leaving the rest not defined
        """
        self.BCTypeCheck(BCType)
        if ("Advective" in BCType) and ("vof" not in BCType):
            BC=self.constantValue(0.0)
        else:
            return self.empty()



    def noSlip(self,BCType):
        """Imposes a noSlipCondition
        Takes the following arguments
        BCType: Type of boundary condition 
        Returns Dirichlet and Advective conditions as zero, leaving the rest not defined
        """
        self.BCTypeCheck(self,BCType)
        if (BCType is "uDirichlet") or (BCType is "vDirichlet") or (BCType is "wDirichlet") or ("Advective" in BCType):
            BC=self.constantValue(0.0)
            if "vof" in BCType:
                BC = self.empty()
            return BC
        else:
            return self.empty()
#: U=0
    def twoPhaseVelocityInlet(self,BCType,x,U,seaLevel,b_or,vert_axis=1,air=1.0,water=0.0):
        """Imposes a velocity profile at the lower phase and an open boundary at the upper face
        Takes the following arguments
        BCType: Type of boundary condition 
        x: Position vector
        U: Velocity vector at the global system
        seaLevel: water level at global coordinate system
        b_or is the orientation of the boundary (always positive outwards vector)
        vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1]
        air: Volume fraction for air (1.0 by default)
        water: Volume fraction for water (
        Below the seawater level, the condition returns the Dirichlet and pAdvective condition according to the inflow velocity
        Above the sea water level, the condition returns the gravity as zero, and sets Dirichlet condition to zero, only if there is an 
        zero inflow velocity component
        THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
        """
        self.BCTypeCheck(BCType)
        if len(U)!=3: 
            print("Velocity in twoPhaseInflow Boundary condition must have three components")
            exit(1)
        u = float(U[0])
        v = float(U[1])
        w = float(U[2])
# This is the normal velocity, based on the inwards boundary orientation -b_or
        u_p =u*b_or[0]+v*b_or[1]+w*b_or[2]
        u_p = -u_p
        if x[vert_axis]<seaLevel:
            if BCType is "uDirichlet":
                return self.constantValue(u)
            elif BCType is "vDirichlet":
                return self.constantValue(v)
            elif BCType is "wDirichlet":
                return self.constantValue(w)
            elif BCType is "pAdvective":
                return self.constantValue(u_p)
            elif BCType is "vofDirichlet":
                return self.constantValue(water) 
            else: return self.empty()
  
        elif x[vert_axis]>=seaLevel:
            if (BCType is "uDirichlet") and( u==0.):
                return self.constantValue(0.0)
            elif (BCType is "vDirichlet") and (v==0.):
                return self.constantValue(0.0)
            elif (BCType is "wDirichlet") and (w==0.):
                return self.constantValue(0.0)
            elif "Diffusive" in BCType:
                return self.constantValue(0.0)
            elif BCType is "vofDirichlet":
                return self.constantValue(air) 
            else: return self.empty()
            

    def hydrostaticPressureOutlet(self,BCType,rho,g,refLevel,b_or,pRef=0.0,vert_axis=1,air = 1.0):
        """Imposes a hydrostatic pressure profile and  open boundary conditions
        Takes the following arguments
        BCType: Type of boundary condition 
        x: Position vector
        rho: Phase density
        g: Gravitational acceleration vector
        refLevel: Level at which pressure = pRef
        pRef: Reference value for the pressure at x[vert_axis]=refLevel, be default set to 0
        b_or is the orientation of the boundary (always positive outwards vector)
        vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1
        Returns a hydrostatic profile for pDirichlet and zero Diffusive conditions. 
        If the boundary is aligned with one of the main axes, sets the tangential velocity components to zero as well
        THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
        """
        self.BCTypeCheck(BCType)
        a0 = pRef - rho*g[vert_axis]*refLevel
        a1 = rho*g[vert_axis]
# This is the normal velocity, based on the boundary orientation
        if BCType is "pDirichlet":
            return self.linear(a0,a1,vert_axis)
        elif (BCType is "uDirichlet") and( b_or[0]==0):
            return self.constantValue(0.)
        elif (BCType is "vDirichlet") and (b_or[1]==0):
            return self.constantValue(0.)
        elif (BCType is "wDirichlet") and (b_or[2]==0):
            return self.constantValue(0.)
        elif "Diffusive" in BCType:
            return self.constantValue(0.0)        
        elif BCType is "vofDirichlet":
            return self.constantValue(air)        
        else: return self.empty()
               
# Derived boundary conditions
    def hydrostaticPressureOutletWithDepth(self,BCType,x,seaLevel,rhoUp,rhoDown,g,refLevel,b_or,pRef=0.0,vert_axis=1,air=1.0,water=0.0):
        """Imposes a hydrostatic pressure profile and open boundary conditions with a known otuflow depth
        Takes the following arguments
        BCType: Type of boundary condition 
        x: Position vector
        rhoUp: Phase density of the upper part
        rhoDown: Phase density of the lower part
        g: Gravitational acceleration vector
        refLevel: Level at which pressure = pRef
        pRef: Reference value for the pressure at x[vert_axis]=refLevel, be default set to 0
        b_or is the orientation of the boundary (always positive outwards vector)
        vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1
        Returns the hydrostaticPressureOutlet except when the pressure and the vof are defined. Then it returns the pressure 
        and vof profile based on the known depth
        If the boundary is aligned with one of the main axes, sets the tangential velocity components to zero as well
        THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
        """
        self.BCTypeCheck(BCType)
        BC = self.hydrostaticPressureOutlet(BCType,rhoUp,g,refLevel,b_or,pRef,vert_axis,air)
        if BCType is "pDirichlet" and x[vert_axis]<seaLevel:
            a0 = pRef - rhoUp*g[vert_axis]*(refLevel - seaLevel) - rhoDown*g[vert_axis]*seaLevel
            a1 = rhoDown*g[vert_axis]
            return self.linear(a0,a1,vert_axis)
        elif BCType is "vofDirichlet" and x[vert_axis]<seaLevel:
            return self.constantValue(water)
        else: return BC
     def forcedOutlet(self,BCType,x,U,seaLevel,rhoUp,rhoDown,g,refLevel,b_or,pRef=0.0,vert_axis=1,air=1.,water=0.):
        """Imposes a known velocit & pressure profile at the outflow
        Takes the following arguments
        BCType: Type of boundary condition 
        BCType: Type of boundary condition 
        x: Position vector
        rhoUp: Phase density of the upper part
        rhoDown: Phase density of the lower part
        g: Gravitational acceleration vector
        refLevel: Level at which pressure = pRef
        pRef: Reference value for the pressure at x[vert_axis]=refLevel, be default set to 0
        b_or is the orientation of the boundary (always positive outwards vector)
        vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1
        Returns only the dirichlet condition usinghydrostaticPressureOutletWithDepth for pressure and twoPhaseVelocityInlet 
        for velocity, wuth inversed velocity
        THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
        """
        self.BCTypeCheck(BCType)
        BC1 = self.hydrostaticPressureOutletWithDepth(BCType,x,seaLevel,rhoUp,rhoDown,g,refLevel,b_or,pRef,vert_axis,air,water)
        BC2 = self.twoPhaseVelocityInlet(BCType,x,U,seaLevel,b_or,vert_axis,air,water)
        if "Dirichlet" in BCType:
            if BCType is ("uDirichlet" or "vDirichlet" or "wDirichlet" or "vofDirichlet"):
                return BC2
            else: return BC1

     def regularWaves(self)
        """This condition is to generate a regular wave at an inlet
           1st order theory is used by default and FFT wave theory is used if an array of coefficients is provided
           The conditions is derived from the velocity inlet 
        """
        z = x[vert_axis] - seaLevel
        self.BCTypeCheck(BCType)
        k = 2*pi/lw
        kdir = k*wdirection
        sigma = 2*pi/T
        phi = sum(kdir * x) - omega*t + phi0
        Ncoeff = len(Y)
        if(wdirection[vert_axis]!=0):
            print("Gravity axis is aligned with a directional wave component. Check input")
            exit(1)    
        if(waveType == "linear"):
            B[0] = 0.5*height*ww
            Y[0] = 0.5*height
            B[1:] = 0.
            Y[1:] = 0.
        elif(waveType == "nonlinear"):
            for i in range(Y):
                Y[i]/=k
            for i in range(B):
                B[i]*=sqrt(abs(g[vert_axis])/k)*(i+1)
        for 
        
def waveVelocity_u(x,t):
     return sigma*amplitude*cosh(k*(z(x)+h))*cos(theta(x,t))/sinh(k*h)

def waveVelocity_v(x,t):
     return sigma*amplitude*sinh(k*(z(x)+h))*sin(theta(x,t))/sinh(k*h)

 
#solution variables

def wavePhi(x,t):
    return x[1] - waveHeight(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_v(x,t):
    waterspeed = waveVelocity_v(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[1]+(1.0-H)*waterspeed

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

outflowHeight=inflowHeightMean

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[1] - outflowHeight)

def outflowPhi(x,t):
    return x[1] - outflowHeight

def outflowPressure(x,t):
  if x[1]>inflowHeightMean:
    return (L[1]-x[1])*rho_1*abs(g[1])
  else:
    return (L[1]-inflowHeightMean)*rho_1*abs(g[1])+(inflowHeightMean-x[1])*rho_0*abs(g[1])


    #p_L = L[1]*rho_1*g[1]
    #phi_L = L[1] - outflowHeight
    #phi = x[1] - outflowHeight
    #return p_L -g[1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
    #                                                     -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

def twpflowVelocity_w(x,t):
    return 0.0

def zeroVel(x,t):
    return 0.0

            

                

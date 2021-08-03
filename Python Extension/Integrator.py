class Integrator:
    
    import math
    import numpy as np
    import Particle as p


    #initializing shared variables
    q0 = None
    q1 = None
    p1 = None
    f1 = None
    ee = None
    h2 = None
    h12 = None
    eps = math.pow(10, -6)
    '''
    def __init__(self):
        pass
    '''
    
    def __init__(self, xmax, nx, e1, delta_E):
        self.xmax = xmax
        self.nx = nx 
        self .e1 = e1
        self.delta_E = delta_E
    
    

    def setEpsilon(self, eps):
        self.eps = eps 

    def setRange(self, xmax):
        self.xmax = xmax

    def setNumOfPoints(self, nx):
        self.nx = nx
    
    def setInitialEnergy(self, e1):
        self.e1 =e1
    
    def setDeltaE(self, delta_E):
        self.delta_E = delta_E

    def potential(self, x):
        #infinite square well potential
        return 0.0

        #square well with walls V(|x| >1) = 10
        #if abs(x) < 1:
        #    return 0.0
        #else:
        #    return 10.0
    
        #infinite square well with a box between -0.5 and 0.5
        #   return 5.0
        #else:
        #    return 0.0

        #simple harmonic oscillator potential


    def normalize(self, h):
        

        norm = self.math.pow(self.psi[0], 2) + self.math.pow(self.psi[-1], 2)

        for i in range(1, len(self.psi)-3, 2):

            norm = norm + 4.0 * self.math.pow(self.psi[i], 2) +2.0 * self.math.pow(self.psi[i+1], 2)

        norm = norm + 4.0 * self.math.pow(self.psi[-2], 2)

        norm = 1.0/self.math.sqrt(norm*h/3.0)

        for i in range(len(self.psi)):

            self.psi[i] = self.psi[i] * norm


    #integrates one step using Numerov's method. 
    def numerovstep(self, x):
        

        q2 = self.h2 * self.f1 * self.p1 + 2.0 * self.q1 - self.q0

        self.q0 = self.q1; self.q1 = q2

        self.f1 = 2.0*(self.potential(x)-self.ee)

        self.p1 = self.q1/(1.0 - self.h12 * self.f1)


    def setinitcond(self,h): #sets the initial boundary conditions


        self.psi[0] = 0.0
        self.psi[1] = 0.00010

        self.p1 = self.psi[0]
        self.f1 = 2.0  * (self.potential(-self.xmax)-self.ee)        
        self.q0 = self.psi[0] * (1.0 - self.h12 * self.f1)
        self.f1 = 2.0  * (self.potential(-self.xmax+h)-self.ee)        
        self.q1 = self.psi[1] * (1.0 - self.h12 * self.f1)

       


    def integrate(self, e1): #runs the integration
        


        self.ee = e1

        h = self.xmax/self.nx

        self.h2 = self.math.pow(h, 2)

        self.h12 = self.h2/12.0

        self.setinitcond(h) 

        for i in range(2, int(len(self.psi)/2)):
            
            x = (i-self.nx)  * h

            self.numerovstep(x)

            self.psi[i] = self.p1

            self.psi[len(self.psi)-i] = self.psi[i]

        

        
        self.normalize(h)
        
    
    def writemessage(self, n, e, b):

        print(n,"  E =", e, "Boundary deviation " , b)
        
    
    
    def findWavefunction(self):

        
        #variable initializatio
   

        '''
        xmax = float(input("Enter the maximum value of x: ")) 
        nx = int(input("Enter the number of points either side of 0: "))
        e1 = float(input("Searching of E>E0, give E0: "))
        delta_E = float(input("Delta-E in initial course search: "))
        eps = float(input("enter boundary convergence parameter"))

        '''

        self.psi = self.np.empty((2 * self.nx)+1)
    

        #starting course search
        ni = 0
        #print("starting course search")
        self.integrate(self.e1)
        b1 = self.psi[-1]
        #self.writemessage(ni, self.e1, b1)

    
        while True:
            ni += 1

            e2 = self.e1 + self.delta_E

            self.integrate(e2)

            b2 = self.psi[-1]

            #self.writemessage(ni, e2, b2)
        

            if(b1 * b2 < 0.0 ): break


            self.e1 = e2
            b1 = b2

        ni = 0
    
        #print("starting bisection")
    
        
        while True:

            ni += 1
            e0 = (self.e1 + e2)/2.0

            self.integrate(e0)

            b0 = self.psi[-1]
        

            #self.writemessage(ni, e0, b0)

            if(abs(b0) <= self.eps): break

        

            if(b0 * b1 <= 0.0):
            
                e2 = e0
                b2 = b0

            else:

                self.e1 = e0
                b1 = b0



        '''
        output = open("Psi.dat", "w")


        for i in range(len(self.psi)):
        
            output.write(str((i-self.nx)*self.xmax/self.nx)+'      '+ str(self.psi[i])+"\n")

        output.close()
        '''

        output = self.p.Particle(e0, self.psi)

        return output


'''
For SHO and square well, find even parity solutions find the center point 
at x=0 and then check if the derivative is zero. then reflect the 
function around the origin. You can check the derivative around 0 using 
the difference method.

For Odd parity solutions check psi itself it should be 0 at x = 0

After changing the potential, pretty much everything else can stay the 
same, setinitcond will work for a SHO with sufficiently large xmax
'''

#solving the 1D Schrodinger Equation using the Shooting method and Numerov's algorithm

import math
import numpy as np

global eps
global q0
global q1
global p1
global f1
global ee
global h2
global h12

eps = math.pow(10, -6)


def potential(x): #calculates the potential as a function of x. 

    
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


def normalize(nx, numOfPoints, h, psi):
    
    global eps
    global q0
    global q1
    global p1
    global f1
    global ee
    global h2
    global h12
    

    norm = math.pow(psi[0], 2) + math.pow(psi[numOfPoints-1], 2)

    for i in range(1, len(psi)-3, 2):

        norm = norm + 4.0 * math.pow(psi[i], 2) +2.0 * math.pow(psi[i+1], 2)

    norm = norm + 4.0 * math.pow(psi[numOfPoints-2], 2)

    norm = 1.0/math.sqrt(norm*h/3.0)

    for i in range(len(psi)):

        psi[i] = psi[i] * norm


#integrates one step using Numerov's method. 
def numerovstep(x):
    
    global eps
    global q0
    global q1
    global p1
    global f1
    global ee
    global h2
    global h12

    q2 = h2 * f1 * p1 + 2.0 * q1 - q0

    q0 = q1; q1 = q2

    f1 = 2.0*(potential(x)-ee)

    p1 = q1/(1.0 -h12 * f1)


def setinitcond(xmax, h, psi0, psi1): #sets the initial boundary conditions

    global eps
    global q0
    global q1
    global p1
    global f1
    global ee
    global h2
    global h12

    psi0 = 0.0
    psi1 = 0.00010

    p1 = psi0
    f1 = 2.0  * (potential(-xmax)-ee)        
    q0 = psi0 * (1.0 - h12 * f1)
    f1 = 2.0  * (potential(-xmax+h)-ee)        
    q1 = psi1 * (1.0 - h12 * f1)


def integrate(xmax, nx, numOfPoints, e1, psi): #runs the integration
    
    global eps
    global q0
    global q1
    global p1
    global f1
    global ee
    global h2
    global h12


    ee = e1

    h = xmax/nx

    h2 = math.pow(h, 2)

    h12 = h2/12.0

    setinitcond(xmax, h, psi[0], psi[1]) 

    for i in range(2, len(psi)):
        
        x = i  * h

        numerovstep(x)

        psi[i] = p1

        #print(i)

    
    normalize(nx,numOfPoints, h, psi)
   


def writemessage(n, e, b):
    
    print(n,"  E =", e, "Boundary deviation " , b)

 
def main():
    
   #variable initialization
    global eps
    global q0
    global q1
    global p1
    global f1
    global ee
    global h2
    global h12
    e0 = 0.0
    e1 = 0.0
    e2 = 0.0
    b0 = 0.0
    b2 = 0.0

    xmax = float(input("Enter the maximum value of x: ")) 
    nx = int(input("Enter the numbe of points either side of 0: "))
    e1 = float(input("Searching of E>E0, give E0: "))
    delta_E = float(input("Delta-E in initial course search: "))


    numOfPoints = 2 * nx

    psi = np.empty(numOfPoints)
    

    #starting course search
    ni = 0
    print("starting course search")
    integrate(xmax, nx, numOfPoints, e1, psi)
    b1 = psi[numOfPoints-1]
    writemessage(ni, e1, b1)


    while True:
        ni += 1

        e2 = e1 + delta_E

        integrate(xmax, nx, numOfPoints, e2, psi)

        b2 = psi[numOfPoints-1]

        writemessage(ni, e2, b2)
        #print(b1)

        if(b1 * b2 < 0.0 ): break


        e1 = e2
        b1 = b2

    ni = 0

    print("starting bisection")

    while True:

        ni += 1
        e0 = (e1 + e2)/2.0

        integrate(xmax, nx, numOfPoints, e2, psi)

        b0 = psi[numOfPoints-1]

        writemessage(ni, e0, b0)

        if(abs(b0) <= eps): break

        if(b0 * b1 <= 0.0):
            
            e2 = e0
            b2 = b0

        else:

            e1 = e0
            b1 = b0

    #j =0
    #while True:
    #    print(j)
    #    print(psi[j])
    #    j+=1

    output = open("myPsi.dat", "w")
    print(len(psi))

    for i in range(len(psi)):
        
        output.write(str((i-nx)*xmax/nx)+'      '+ str(psi[i])+"\n")

    output.close()
            

if __name__ == "__main__":
    main()






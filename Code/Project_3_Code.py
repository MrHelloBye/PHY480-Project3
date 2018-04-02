# PHY480 - Project 3
# Planet Class and Integrators
# Daniel Coulter, David Butts

import numpy as np

class Planet:
    """
    A generalized planet object that has a circular orbit, location, orbital
    period, and mass.
    """
    
    def __init__(self, mass, xloc, yloc, period, sun=False):
        """
        Creates a planet object with a mass and location in the x-y plane.
        :param mass: Mass of the planet in solar masses
        :param xloc: x location in the x-y plane in AU
        :param yloc: y location in the x-y plane in AU
        """
        self.mass = mass
        self.xloc = xloc
        self.yloc = yloc        
        if sun == True:
            self.xvel = 0.0
            self.yvel = 0.0
        else:
            self.xvel = 0.0
            self.yvel = 2*np.pi*xloc/period
        
        
    def location(self):
        """
        Returns the current coordinates of a planet for reference.
        """
        
        coordinates = [self.xloc, self.yloc]
              
        return coordinates
    

def force(p1,p2):
    """
    Calculates the force between two objects.
    :param p1: First of two planet objects
    :param p2: Second of two planet objects
    :return: x and y accelerations of p1 as a list
    """
    
    xsep = p1.xloc-p2.xloc
    ysep = p1.yloc-p2.yloc
    r = np.sqrt(xsep**2+ysep**2)
    tforce = (4*np.pi**2)*p1.mass*p2.mass/(r**2)
    
    if xsep == 0:
        xforce = 0
    else:
        xforce = -tforce*xsep/r
    
    if ysep == 0:
        yforce = 0
    else:
        yforce = -tforce*ysep/r
    
    xaccel = xforce/p1.mass
    yaccel = yforce/p1.mass
    
    return [xaccel,yaccel]


def forward_euler(p1, p2, h):
    """
    Solves differential equations using the forward Euler method.
    :param p1: First of two planet objects
    :param p2: Second of two planet objects
    :param h: Step size
    """
    
    p1.xvel += h*force(p1,p2)[0]
    p1.yvel += h*force(p1,p2)[1]
    
    p2.xvel -= h*force(p2,p1)[0]
    p2.yvel -= h*force(p2,p1)[1]
    
    p1.xloc += h*p1.xvel
    p1.yloc += h*p1.yvel
    
    p2.xloc -= h*p2.xvel
    p2.yloc -= h*p2.yvel
    

def velocity_verlet(pl, h):
    """
    Solves differential equations using Velocity Verlet method.
    :param pl: A list of the planets to be used
    :param h: Step size
    """
    
    xpos = []
    ypos = []
    a_x = []
    a_y = []
    new_a_x = []
    new_a_y = []
    for p in pl:
        xpos.append(p.xloc)
        ypos.append(p.yloc)
        a_x.append(0)
        a_y.append(0)
        new_a_x.append(0)
        new_a_y.append(0)
        
    for a in range(0, len(pl)):
        for b in range(0, len(pl)):
            if a != b:
                xacc = force(pl[a],pl[b])[0]
                yacc = force(pl[a],pl[b])[1]
                a_x[a] += xacc
                a_y[a] += yacc
    
    for a in range(0, len(pl)):
        pl[a].xloc += h*pl[a].xvel + (h**2/2)*a_x[a]
        pl[a].yloc += h*pl[a].yvel + (h**2/2)*a_y[a]
        
    for a in range(0, len(pl)):
        for b in range(0, len(pl)):
            if a != b:
                new_xacc = force(pl[a],pl[b])[0]
                new_yacc = force(pl[a],pl[b])[1]
                new_a_x[a] += new_xacc
                new_a_y[a] += new_yacc
                
    for a in range(0, len(pl)):
        pl[a].xvel += (h/2) * (a_x[a]+new_a_x[a])
        pl[a].yvel += (h/2) * (a_y[a]+new_a_y[a]) 


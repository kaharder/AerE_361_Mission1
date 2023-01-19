# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:30:23 2020

@author: Keith Harder
"""

#For my modifcations I:
#made the camera track the rocket instead of a fixed view.
#changed the rocket specs to those of a falcon heavy.
#made the rocket two stage, with the line shortening when the first stage is dropped.
#made the time step variable. when 'm' is pressed, the timestep increases.
#when 'n' is pressed, the timestep decreases. Negative time is weird but fun.
#fixed the proglem where we were never told to make the distance between the earth and rocket change.
#made the rocket "dissapear" when it hits the earth; computations still run but arrow is infinismal. 
#cleaned up code, combined individual if statements

# Template begins here
# --------------------

# Include import statements and unit registry from above

# (Remove the above line for Spyder; 
# you will need to type '%matplotlib qt' into the command window instead.)
import numpy as np
from matplotlib import pyplot as pl
import pint
ur = pint.UnitRegistry() # Import pint and create unit registry

# Include Rocket and planet characteristics; physics constants

mempty = 72300*ur.kg
mp = 16800*ur.kg
mf = 1340200*ur.kg
mtot = mp+mf+mempty
me = 5.98e24*ur.kg
tsfc = 0.00038*(ur.s/ur.m)
thrust = 68457000*ur.N 
G = 6.674e-11*(ur.m**3 / (ur.kg * ur.s**2))
dt = 0.0*ur.s

# Include earth diameter, surface speed, initial rocket position,
# initial rocket velocity

# Earth's diameter is given by
d_earth = 12.7e6*ur.m 
# and Earth's tangential speed (due to rotation) is 
# approximately
v_surface = 460.0*ur.m/ur.s # at the equator

# Define initial thrust angle and direction

# Position vector
# This puts the horizontal (x) position at zero
# and the vertical (y) position at one earth radius
# (will be the initial position of our rocket)
pos = np.array((0.0,d_earth.to(ur.m).magnitude/2.0))*ur.m


#created array vel to store velocity and gave it an initial velocity equal to the speed of the surface of the earth in the x
vel = np.array((v_surface.to(ur.m/ur.s).magnitude, 0.0))*ur.m/ur.s


# Given, for example
T_angle = 90*ur.degrees # 60 deg. CCW from horizontal

# Evaluate and store T_direc. 
# The components of T_direc are the sine and cosine of T_angle
# (but you have to figure out which is which and whether there
# are any minus signs). Also the np.sin() and np.cos() functions
# take their parameters in radians, not degrees. 

# +x should be to the right and +y should be up. Make sure the 
# direction vector makes sense. 

# T_direction = 
T_direc = np.array((np.cos(T_angle.to(ur.radians).magnitude), np.sin(T_angle.to(ur.radians).magnitude)))


# Include code for the baseline plot and initial arrow

pl.ion() # interactive mode ON
fig = pl.figure()
pl.axis((-1e6,1e6,6e6,7e6)) # Define axis bounds in (assumed) meters
pl.grid(True) # Enable grid

# create green filled circle representing earth
earth = pl.Circle((0,0),float(d_earth/2.0/ur.m),color='g')
fig.gca().add_artist(earth)

# Add your Arrow-drawing code here
arrowplt = pl.arrow(pos[0].to(ur.m).magnitude, pos[1].to(ur.m).magnitude, T_direc[0]*300000, T_direc[1]*300000, width= 90000,head_width= 90000, head_length= 150000)
arrowplt
#the scale seemed too big when I ploted the rocket the first time, so I shrunk it. If it ends up being too small later, I can grow it again.

# Include the event handler 

def event_handler(event):
    global dt
    if event.key=='m':
        dt = dt+.2*ur.s
        pass
    elif event.key=='n':
        dt = dt-.2*ur.s
        pass
    global T_angle
    if event.key==',': # press comma to rotate left
        # Rotate left (increase angle)
        T_angle += 10.0*ur.degree
        pass
    elif event.key=='.': # press period to rotate right
        # Rotate right (decrease angle)
        T_angle -= 10.0*ur.degree
        pass
    elif hasattr(event,"button") and event.button==1:
        if event.xdata < 0.0: 
            # Click on left-half of plot to
            # Rotate left (increase angle)
            T_angle += 10.0*ur.degree
            pass
        else:
            # Click on right-half of plot to
            # Rotate right (decrease angle)
            T_angle -= 10.0*ur.degree
            pass
        pass
    
    pass

# Connect this event handler to the figure so it is called
# when the mouse is clicked or keys are pressed. 
fig.canvas.mpl_connect('key_press_event',event_handler)
fig.canvas.mpl_connect('button_press_event',event_handler)


# Include the vector norm function from the 
# force of gravity calculation, above. 

# In the last lab we used np.linalg.norm() to find the magnitude of a vector.
# Unfortunately that function does not work with quantities that have units. 
# Instead, use this function norm() that does work with quantities that have units
def norm(vec):
    """Unit aware vector norm (magnitude)"""
    return np.sqrt(np.sum(vec**2.0))


# r =  # Should have units of distance (meters)
# F_Gravity = # Should have units of Force (N) or equivalent
r = np.sqrt((pos[0]**2)+(pos[1]**2))
print(r)
F_Gravity = (G*(mtot)*(me))/(r**2)
print(F_Gravity)
print(mtot)

t=0.0 * ur.s # Start at t=0 

# Loop until ctrl-C or press stop button on 
# Spyder console or t exceeds 36000 seconds (10 hours
# of simulated time)

while t < 36000 * ur.s: 
    if(r.magnitude<6350000):
       hl = 0
       w=0
       l=0
       pass
    else:
        hl=150000
        w=90000
        l=100000
        pass
    
    arrowplt.remove()
    
    if(mf>(107500*ur.kg)):
        T_Vec = T_direc*thrust
        mf=mf-(tsfc*thrust*dt)
        arrowplt = pl.arrow(pos[0].to(ur.m).magnitude, pos[1].to(ur.m).magnitude, T_direc[0]*3*l, T_direc[1]*3*l, width= w,head_width= w, head_length= hl)
        
    elif(mf> (0*ur.kg)):
        mempty = 4000*ur.kg
        mtot = mp+mf+mempty
        thrust = 981000*ur.N  
        tsfc = 0.000276*(ur.s/ur.m)
        T_Vec = T_direc*thrust
        mf=mf-(tsfc*thrust*dt)
        arrowplt = pl.arrow(pos[0].to(ur.m).magnitude, pos[1].to(ur.m).magnitude, T_direc[0]*l, T_direc[1]*l, width= w,head_width= w, head_length= hl)
        
    elif(mf<(0*ur.kg)):
        mf=0*ur.kilograms
        arrowplt = pl.arrow(pos[0].to(ur.m).magnitude, pos[1].to(ur.m).magnitude, T_direc[0]*l, T_direc[1]*l, width= w,head_width= w, head_length= hl)

    else:
        T_Vec = T_direc*0.0*ur.newton
        mf=0*ur.kilograms
        arrowplt = pl.arrow(pos[0].to(ur.m).magnitude, pos[1].to(ur.m).magnitude, T_direc[0]*l, T_direc[1]*l, width= w,head_width= w, head_length= hl)

        
    pl.axis((pos[0].magnitude-3e6, pos[0].magnitude+3e6, pos[1].magnitude-3e6, pos[1].magnitude+3e6)) 
    # Update the plot:     
    #  * Call arrowplt.remove() method on the old arrow
    #  * Calculate the thrust (rocket) direction from
    #    the rocket angle (code way above)
    #  * Plot the new arrow (code way above)
    #  * Label the fuel state in the plot title, e.g.
    #    pl.title("Fuel remaining %f kg" % ())
    #  * Show the time in the x axis label, e.g.
    #    pl.xlabel("Time = %f s" % ())
    #  * Select the plot region according to fuel state (code above)
    
    T_direc = np.array((np.cos(T_angle.to(ur.radians).magnitude), np.sin(T_angle.to(ur.radians).magnitude)))
    
    arrowplt
    
    pl.title("Fuel remaining %f kg" % (mf.magnitude))
    
    pl.xlabel("Time = %f s" % (t.to(ur.s).magnitude))
    
    
    # These next two lines cause the plot display to refresh
    fig.canvas.draw()
    fig.canvas.flush_events()

    # Determine the forces on the rocket:  
    #  * Calculate the magnitude of the force of gravity
    #  * Calculate the direction of the force of gravity
    #    and gravity vector on the rocket
    #  * Calculate the thrust vector
    
    # 
    # Direc_Gravity =    # Direction of gravity
    r = np.sqrt((pos[0]**2)+(pos[1]**2))
    F_Gravity = (G*(mtot)*(me))/(r**2)
    Direc_Gravity = np.array((-(pos[0].magnitude)/(r.magnitude), -(pos[1].magnitude)/(r.magnitude)))
    # Vec_Gravity =     # Vector gravitational force on rocket
    Vec_Gravity = Direc_Gravity*F_Gravity
        
    # Determine the change in velocity
    
    F_tot = np.array((T_Vec[0].magnitude+Vec_Gravity[0].magnitude, T_Vec[1].magnitude+Vec_Gravity[1].magnitude))*(ur.kilogram*ur.m)/ur.s**2
    acceleration = (F_tot)/mtot
    dv = acceleration*dt
    vel = np.array((vel[0].magnitude+dv[0].magnitude, vel[1].magnitude+dv[1].magnitude))*ur.m/ur.s

    # Determine the change in position
    
    dpos = (vel)*dt
    pos = np.array((pos[0].magnitude+dpos[0].magnitude, pos[1].magnitude+dpos[1].magnitude))*ur.m
    
    
    t += dt # Update the time
    pass  # End of loop block


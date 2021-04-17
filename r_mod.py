import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from ROHSApy import ROHSA
import astropy
import astropy.units as u
import astropy.constants as c
from spectral_cube import SpectralCube
import time
from fractions import Fraction


def plot_mesh(param, RA, ra1, ra2, DEC, dec1, dec2, data, mask, n_gauss):
    """
    Plots mosiac of one parameter as a pcolormesh, alternating masked and unmasked versions of the data.
    
    
    """

    cmap = ''
    if param == 'integral':
        cmap = 'inferno'
    elif param == 'position':
        cmap = 'coolwarm'
    elif param == 'dispersion':
        cmap = 'cubehelix'
    else:
        raise Exception("Not a valid parameter.")
        
    

    plt.figure(figsize=(16, 6.67*n_gauss)) 
    
    for i in range(1, 2*n_gauss+1):
        plt.subplot(n_gauss, 2, i)
        if i%2 == 0:
            plt.pcolormesh(RA[ra1:ra2], DEC[dec1:dec2], mask*data[(i-1)//2], cmap=cmap)
        else:
            plt.pcolormesh(RA[ra1:ra2], DEC[dec1:dec2], data[(i-1)//2], cmap=cmap)
        plt.colorbar()
        
    plt.show()
    

    
def plot_cont(param, RA, ra1, ra2, DEC, dec1, dec2, data, mask, n_gauss):
    """
    Plots mosiac of one parameter as a contourf, alternating masked and unmasked versions of the data.
    
    
    """

    cmap = ''
    if param == 'integral':
        cmap = 'inferno'
    elif param == 'position':
        cmap = 'coolwarm'
    elif param == 'dispersion':
        cmap = 'cubehelix'
    else:
        raise Exception("Not a valid parameter.")
        
    

    plt.figure(figsize=(16, 6.67*n_gauss)) 
    
    for i in range(1, 2*n_gauss+1):
        plt.subplot(n_gauss, 2, i)
        if i%2 == 0:
            plt.contourf(RA[ra1:ra2], DEC[dec1:dec2], mask*data[(i-1)//2], cmap=cmap)
        else:
            plt.contourf(RA[ra1:ra2], DEC[dec1:dec2], data[(i-1)//2], cmap=cmap)
        plt.colorbar()
        
    plt.show()

    
def vis(g, mask, outline, levels):
    """
    Visualize the Velocity map of a component
    """
    
    plt.figure(figsize=(10, 8))
    plt.contour(RA[ra1:ra2], DEC[dec1:dec2], position[g-1], colors="k", linewidths=1, alpha=0.3, levels=levels)
    plt.pcolormesh(RA[ra1:ra2], DEC[dec1:dec2], position[g-1], cmap="coolwarm", alpha=0.3)
    plt.contour(RA[ra1:ra2], DEC[dec1:dec2], mask*position[g-1], colors="k", linewidths=1, levels=levels)
    plt.pcolormesh(RA[ra1:ra2], DEC[dec1:dec2], mask*position[g-1], cmap="coolwarm")
    plt.colorbar()

    plt.contour(RA[ra1:ra2], DEC[dec1:dec2], outline, colors="white", linewidths=2)
    plt.title('Visualizing the velocity map of g'+str(g))
    plt.show()
    
    
def dist(x1, y1, x2, y2):
    """
    find distance between two points in the cartesian plane
    """

    d = np.sqrt(  (x2-x1)**2 + (y2-y1)**2 )
    return d



def fit_g(data, rad_v, x, y, wid_guess = 1e5, std = 0):
    """
    Fit single gaussian to one spectra, specified by pixel
    
    Parameters:
    
    data - 3d array
    rad_v - array
    x, y - int
    
    Returns:
    
    amp - float
    cen - float
    wid - float
    
    """
    
    spec = m.int_pix(data, x, y)
    sig = np.full((spec.shape), std)
    
    init_vals = [np.max(spec), rad_v[np.where(spec == np.max(spec))[0][0]], wid_guess] # for [amp, cen, wid]
    try:
        best_vals, covar = curve_fit(gauss, rad_v, spec, p0=init_vals, sigma=sig)
    except RuntimeError:
        best_vals = [0, 0, 1]
    amp, cen, wid = best_vals
    
    
    return amp, cen, wid



def gauss(x, amp, cen, wid):
    """
    Initialize a gaussian function
    
    """
    if wid == 0:
        wid = 1
    else:
        wid = wid
        
    return amp * np.exp(-(x-cen)**2 / wid)

def quad(x, a, b, c):
    """
    Initialize a quadratic function
    """
    
    return a*(x-b)**2 + c

def TtF(array):
    """
    turns an array of True/False to False/True
    """
    
    for i in range(len(array)):
        if array[i]:
            array[i] = False
        elif not array[i]:
            array[i] = True
        else:
            continue
    return array
    
    
    











def plot_mosaic():
    """
    Plot a mosaic of the gaussian decomposition across a specified pixel range
    
    Parameters:
    
    
    Returns:
    
    
    
    """
    
    
    #Plot mosaic spectra
    pvalues = np.logspace(-1, 0, 2)
    pmin = pvalues[0]
    pmax = pvalues[-1]

    def norm(pval):
         return (pval - pmin) / float(pmax - pmin)

    ny = 10; nx = 10
    center_y = 32; center_x = 30
    cb = "magenta"
    cw = "crimson"
    fig, axs = plt.subplots(ny, nx, sharex=True, sharey=True,figsize=(10.,6.))
    fig.subplots_adjust(hspace=0, wspace=0, left=0, right=1, top=1, bottom=0)
    for i in np.arange(ny):
        for j in np.arange(nx):
            axs[i][j].step(core.v, data[:,center_y+i,center_x+j], color='cornflowerblue', linewidth=2.)
            axs[i][j].plot(core.v, reconstructed_cube[:,center_y+i,center_x+j], linestyle="-", linewidth=2., color="k")
            for k in range(int(params.shape[0]/3)):
                axs[i][j].plot(core.v, core.gauss(np.arange(len(core.v)), 
                                                gaussian[0::3][k][center_y+i,center_x+j], 
                                                gaussian[1::3][k][center_y+i,center_x+j], 
                                                gaussian[2::3][k][center_y+i,center_x+j]), 
                                                linewidth=2., color=plt.cm.inferno(pvalues[k]))
            if j == 0: axs[i][j].set_ylabel(r'T [k]')
            axs[i][j].set_xlabel(r'v [km s$^{-1}$]')
            
            
            
            
def grad_dec(grid, start):
    """
    Follows the path of quickest descent and records the height at each step to draw slope.
    
    Parameters:
    
    grid - scalar field, array - field which upon which is descended
    start - tuple, list - starting point where the "ball" is dropped
    
    Returns:
    
    heights - list - the gradient of the field
    
    """
    
    loc = start
    
    steps = 0
    smallest = -1
    heights = []
    pt = np.zeros(grid.shape)
    
    while smallest < 0:
        

        pt[loc[1], loc[0]] = 10
        mask = pt > 0
        g1 = np.empty(mask.shape)
        g1.fill(np.nan)
        g1[mask == True] = 1
        
        plt.figure(figsize=(5, 5))
        plt.pcolormesh(grid)
        plt.pcolormesh(g1*pt, cmap='Greys')
        plt.show()

        here = grid[loc[1]][loc[0]]
        heights.append(here)
        
        deltas = []
        
        #Making sure within upper bound
        if loc[1] < grid.shape[1]-1:
            up = grid[loc[1]+1][loc[0]]
            d_up = up-here
            deltas.append(d_up)
        else:
            d_up = 0
            
        #Making sure within lower bound
        if loc[1] > 0:
            down = grid[loc[1]-1][loc[0]]
            d_down = down-here
            deltas.append(d_down)
        else:
            d_down = 0

        #Making sure within left most bound
        if loc[0] > 0:
            left = grid[loc[1]][loc[0]-1]
            d_left = left-here
            deltas.append(d_left)
        else:
            d_left = 0

        #Making sure within right most bound
        if loc[0] < grid.shape[0]-1:
            right = grid[loc[1]][loc[0]+1]
            d_right = right-here
            deltas.append(d_right)
        else:
            d_right = 0


        smallest = np.min(deltas)

        
        if smallest<0:
            if d_up == smallest:
                loc[1] = loc[1]+1
                steps+=1
            if d_down == smallest:
                loc[1] = loc[1]-1
                steps+=1
            if d_left == smallest:
                loc[0] = loc[0]-1
                steps+=1
            if d_right == smallest:
                loc[0] = loc[0]+1
                steps+=1
        else:
            print("local min found")

       

        time.sleep(1/2)

        
    return heights



    
def check_where(loc, grid, d):
    """
    Check if the next point is either boundary or nan
    Return the curve at that time
    
    Parameters:
    
    loc - tuple
    grid - 2d array
    d - string - can either be 'up' or 'right'
    
    Returns:
    
    isbound - bool
    
    """
    
    isbound = False
    
    if d == 'up':
        
        
        if loc[1]+1 > grid.shape[1]-1:
            isbound = True
            return isbound
        
        nex = grid[loc[1]+1, loc[0]]
        
        if np.isnan(nex):
            isbound = True
            return isbound
    elif d == 'down':
        
        
        if loc[1]-1 < 0:
            isbound = True
            return isbound
        
        nex = grid[loc[1]-1, loc[0]]
        
        if np.isnan(nex):
            isbound = True
            return isbound
        
    elif d == 'right':     
        
        if loc[0]+1 > grid.shape[0]-1:
            isbound = True
            return isbound
        
        nex = grid[loc[1], loc[0]+1]
        
        if np.isnan(nex):
            isbound = True
            return isbound
        
    else:
        raise Exception('not a valid direction')
    
    return isbound
    
    
    
                
def line_dec(start, grid, slope, c='viridis', w=False):
    """
    choose a y=mx+b line, and follow it to the boundary, keeping track of the height along the way
    
    Parameters:
    
    start - tuple - starting point on the grid
    grid - scalar field, array
    slope - float
    
    Returns:
    
    curve - array
    stepsize - float
    
    """
    
    loc = start
    
    m = Fraction(slope).limit_denominator()
    y = m.numerator
    x = m.denominator
    
    steps = 0
    curve = []
    pt = np.zeros(grid.shape)
    
    here = grid[loc[1], loc[0]]
    pt[loc[1], loc[0]] = 10
    mask = pt > 0
    g1 = np.empty(mask.shape)
    g1.fill(np.nan)
    g1[mask == True] = 1
    
    if w:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111, projection=w)
        im1 = ax.imshow(grid, cmap=c)
        im2 = ax.imshow(g1*pt, cmap='Greys')
        
        ra = ax.coords[0]
        ra.set_major_formatter('hh:mm:ss')
        ra.set_separator(':')
        dec = ax.coords[1]
        ax.set_xlabel('RA (J2000)')
        ax.set_ylabel('DEC (J2000)')
        plt.show()
        
    else:
        plt.figure(figsize=(5, 5))
        plt.pcolormesh(grid)
        plt.pcolormesh(g1*pt, cmap='Greys')
        plt.show()
    
    unit = check_type(here)
            
    curve.append(float(here/unit))
    time.sleep(1/2)
    
    while 1 > 0:
        
        y = m.numerator
        x = m.denominator
        
        is_y_bigger = np.abs(y)>=np.abs(x)  
        if is_y_bigger:                            #the numerator is bigger
            
            if y < 0:
                sign = -1
            else:
                sign = 1
            
            ratio = np.round(y/x)                  #find step size for y
            y_step = int(ratio)
            x_step = 1

            iterations = np.ceil(y/ratio)          #find how many times ratio divides into y

            
            
            if iterations < x:                     #if it's less than total x steps, go x dir first

                loc = [loc[0]+1, loc[1]]           #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.show()

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                x=x-x_step

                if check_where(loc, grid, 'right'):
                    return curve

            else:                                  #otherwise, continue as normal

                x=x
                
            y_count = np.abs(y)
            x_count = np.abs(x)

            while y_count>0 and x_count>0:                 #we'll stop when we've gone the total distance
                for i in range(sign*y_step):        #do y steps up, show each one
                    loc = [loc[0] , loc[1]+sign*1]
                    here = grid[loc[1], loc[0]]

                    pt[loc[1], loc[0]] = 10
                    mask = pt > 0
                    g1 = np.empty(mask.shape)
                    g1.fill(np.nan)
                    g1[mask == True] = 1

                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.show()
                    unit = check_type(here)
            
                    curve.append(float(here/unit))
                    time.sleep(1/2)
                    
                    if y > 0:

                        if check_where(loc, grid, 'up'):
                            return curve
                    else:
                        
                        if check_where(loc, grid, 'down'):
                            return curve

                loc = [loc[0]+1, loc[1]]      #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.show()
                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y_count=y_count-y_step                     #update total y left to go
                x_count=x_count-x_step                     #update total x left to go
                if y_count < y_step:                 #if total remaining y is less than a total step, then the step becomes the remainder
                    y_step = y_count

                if check_where(loc, grid, 'right'):
                    return curve


        else:                        #the denominator is bigger
            ratio = np.round(x/y)                  #find step size for x
            y_step = 1
            x_step = int(ratio)

            iterations = np.ceil(x/ratio)          #find how many times ratio divides into x
            
            if x < 0:
                sign = -1
            else:
                sign = 1

            if iterations < y:                     #if it's less than total y steps, go y dir first

                loc = [loc[0], loc[1]+1]           #do 1 step up, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.show()

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y=y-y_step

                if check_where(loc, grid, 'up'):
                    return curve

            else:                                  #otherwise, continue as normal

                y=y
                
            y_count = np.abs(y)
            x_count = np.abs(x)

            while y_count>0 and x_count>0:                 #we'll stop when we've gone the total distance
                for i in range(sign*x_step):        #do x steps right, show each one
                    loc = [loc[0]+1 , loc[1]]
                    here = grid[loc[1], loc[0]]

                    pt[loc[1], loc[0]] = 10
                    mask = pt > 0
                    g1 = np.empty(mask.shape)
                    g1.fill(np.nan)
                    g1[mask == True] = 1

                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.show()
                    unit = check_type(here)
            
                    curve.append(float(here/unit))
                    time.sleep(1/2)

                    if check_where(loc, grid, 'right'):
                        return curve

                loc = [loc[0], loc[1]+sign*1]      #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.show()
                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y_count=y_count-y_step                     #update total y left to go
                x_count=x_count-x_step                     #update total x left to go
                if x_count < x_step:                 #if total remaining x is less than a total step, then the step becomes the remainder
                    x_step = x_count

                if check_where(loc, grid, 'up'):
                    return curve
   
    
    
    return curve

def line_dec_inst(start, grid, slope, c='viridis', w=False):
    """
    choose a y=mx+b line, and follow it to the boundary, keeping track of the height along the way
    Only show final graph
    
    Parameters:
    
    start - tuple - starting point on the grid
    grid - scalar field, array
    slope - float
    
    Returns:
    
    curve - array
    stepsize - float
    
    """
    
    loc = start
    
    m = Fraction(slope).limit_denominator()
    y = m.numerator
    x = m.denominator
    
    steps = 0
    curve = []
    pt = np.zeros(grid.shape)
    
    here = grid[loc[1], loc[0]]
    pt[loc[1], loc[0]] = 10
    mask = pt > 0
    g1 = np.empty(mask.shape)
    g1.fill(np.nan)
    g1[mask == True] = 1

    if w:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111, projection=w)
        im1 = ax.imshow(grid, cmap=c)
        im2 = ax.imshow(g1*pt, cmap='Greys')
        
        ra = ax.coords[0]
        ra.set_major_formatter('hh:mm:ss')
        ra.set_separator(':')
        dec = ax.coords[1]
        ax.set_xlabel('RA (J2000)')
        ax.set_ylabel('DEC (J2000)')
        plt.show()
        
    else:
        plt.figure(figsize=(5, 5))
        plt.pcolormesh(grid)
        plt.pcolormesh(g1*pt, cmap='Greys')
        plt.show()
    
    unit = check_type(here)
            
    curve.append(float(here/unit))
    time.sleep(1/2)
    
    while 1 > 0:
        
        y = m.numerator
        x = m.denominator
        
        is_y_bigger = np.abs(y)>=np.abs(x)  
        if is_y_bigger:                            #the numerator is bigger
            
            if y < 0:
                sign = -1
            else:
                sign = 1
            
            ratio = np.round(y/x)                  #find step size for y
            y_step = int(ratio)
            x_step = 1

            iterations = np.ceil(y/ratio)          #find how many times ratio divides into y

            
            
            if iterations < x:                     #if it's less than total x steps, go x dir first

                loc = [loc[0]+1, loc[1]]           #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                x=x-x_step

                if check_where(loc, grid, 'right'):
                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.show()
                    return curve

            else:                                  #otherwise, continue as normal

                x=x
                
            y_count = np.abs(y)
            x_count = np.abs(x)

            while y_count>0 and x_count>0:                 #we'll stop when we've gone the total distance
                for i in range(sign*y_step):        #do y steps up, show each one
                    loc = [loc[0] , loc[1]+sign*1]
                    here = grid[loc[1], loc[0]]

                    pt[loc[1], loc[0]] = 10
                    mask = pt > 0
                    g1 = np.empty(mask.shape)
                    g1.fill(np.nan)
                    g1[mask == True] = 1
                    
                    unit = check_type(here)
            
                    curve.append(float(here/unit))
                    time.sleep(1/2)
                    
                    if y > 0:

                        if check_where(loc, grid, 'up'):
                            if w:
                                fig = plt.figure(figsize=(5,5))
                                ax = fig.add_subplot(111, projection=w)
                                im1 = ax.imshow(grid, cmap=c)
                                im2 = ax.imshow(g1*pt, cmap='Greys')

                                ra = ax.coords[0]
                                ra.set_major_formatter('hh:mm:ss')
                                ra.set_separator(':')
                                dec = ax.coords[1]
                                ax.set_xlabel('RA (J2000)')
                                ax.set_ylabel('DEC (J2000)')
                                plt.show()

                            else:
                                plt.figure(figsize=(5, 5))
                                plt.pcolormesh(grid)
                                plt.pcolormesh(g1*pt, cmap='Greys')
                                plt.show()
                            return curve
                    else:
                        
                        if check_where(loc, grid, 'down'):
                            if w:
                                fig = plt.figure(figsize=(5,5))
                                ax = fig.add_subplot(111, projection=w)
                                im1 = ax.imshow(grid, cmap=c)
                                im2 = ax.imshow(g1*pt, cmap='Greys')

                                ra = ax.coords[0]
                                ra.set_major_formatter('hh:mm:ss')
                                ra.set_separator(':')
                                dec = ax.coords[1]
                                ax.set_xlabel('RA (J2000)')
                                ax.set_ylabel('DEC (J2000)')
                                plt.show()

                            else:
                                plt.figure(figsize=(5, 5))
                                plt.pcolormesh(grid)
                                plt.pcolormesh(g1*pt, cmap='Greys')
                                plt.show()
                            return curve

                loc = [loc[0]+1, loc[1]]      #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y_count=y_count-y_step                     #update total y left to go
                x_count=x_count-x_step                     #update total x left to go
                if y_count < y_step:                 #if total remaining y is less than a total step, then the step becomes the remainder
                    y_step = y_count

                if check_where(loc, grid, 'right'):
                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.show()
                    return curve


        else:                        #the denominator is bigger
            ratio = np.round(x/y)                  #find step size for x
            y_step = 1
            x_step = int(ratio)

            iterations = np.ceil(x/ratio)          #find how many times ratio divides into x
            
            if x < 0:
                sign = -1
            else:
                sign = 1

            if iterations < y:                     #if it's less than total y steps, go y dir first

                loc = [loc[0], loc[1]+1]           #do 1 step up, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y=y-y_step

                if check_where(loc, grid, 'up'):
                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.show()
                    return curve

            else:                                  #otherwise, continue as normal

                y=y
                
            y_count = np.abs(y)
            x_count = np.abs(x)

            while y_count>0 and x_count>0:                 #we'll stop when we've gone the total distance
                for i in range(sign*x_step):        #do x steps right, show each one
                    loc = [loc[0]+1 , loc[1]]
                    here = grid[loc[1], loc[0]]

                    pt[loc[1], loc[0]] = 10
                    mask = pt > 0
                    g1 = np.empty(mask.shape)
                    g1.fill(np.nan)
                    g1[mask == True] = 1

                    unit = check_type(here)
            
                    curve.append(float(here/unit))
                    time.sleep(1/2)

                    if check_where(loc, grid, 'right'):
                        if w:
                            fig = plt.figure(figsize=(5,5))
                            ax = fig.add_subplot(111, projection=w)
                            im1 = ax.imshow(grid, cmap=c)
                            im2 = ax.imshow(g1*pt, cmap='Greys')

                            ra = ax.coords[0]
                            ra.set_major_formatter('hh:mm:ss')
                            ra.set_separator(':')
                            dec = ax.coords[1]
                            ax.set_xlabel('RA (J2000)')
                            ax.set_ylabel('DEC (J2000)')
                            plt.show()

                        else:
                            plt.figure(figsize=(5, 5))
                            plt.pcolormesh(grid)
                            plt.pcolormesh(g1*pt, cmap='Greys')
                            plt.show()
                        return curve

                loc = [loc[0], loc[1]+sign*1]      #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y_count=y_count-y_step                     #update total y left to go
                x_count=x_count-x_step                     #update total x left to go
                if x_count < x_step:                 #if total remaining x is less than a total step, then the step becomes the remainder
                    x_step = x_count

                if check_where(loc, grid, 'up'):
                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.show()
                    return curve
   
    
    
    return curve

def line_dec_save(start, grid, slope, c='viridis', w=False):
    """
    choose a y=mx+b line, and follow it to the boundary, keeping track of the height along the way, saving each image
    
    Parameters:
    
    start - tuple - starting point on the grid
    grid - scalar field, array
    slope - float
    
    Returns:
    
    curve - array
    stepsize - float
    
    """
    count=0
    
    loc = start
    
    m = Fraction(slope).limit_denominator()
    y = m.numerator
    x = m.denominator
    
    steps = 0
    curve = []
    pt = np.zeros(grid.shape)
    
    here = grid[loc[1], loc[0]]
    pt[loc[1], loc[0]] = 10
    mask = pt > 0
    g1 = np.empty(mask.shape)
    g1.fill(np.nan)
    g1[mask == True] = 1

    if w:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111, projection=w)
        im1 = ax.imshow(grid, cmap=c)
        im2 = ax.imshow(g1*pt, cmap='Greys')
        
        ra = ax.coords[0]
        ra.set_major_formatter('hh:mm:ss')
        ra.set_separator(':')
        dec = ax.coords[1]
        ax.set_xlabel('RA (J2000)')
        ax.set_ylabel('DEC (J2000)')
        plt.savefig('dec_'+str(count)+'.png')
        count+=1
        plt.show()
        
    else:
        plt.figure(figsize=(5, 5))
        plt.pcolormesh(grid)
        plt.pcolormesh(g1*pt, cmap='Greys')
        plt.savefig('dec_'+str(count)+'.png')
        count+=1
        plt.show()
    
    
    unit = check_type(here)
            
    curve.append(float(here/unit))
    time.sleep(1/2)
    
    while 1 > 0:
        
        y = m.numerator
        x = m.denominator
        
        is_y_bigger = np.abs(y)>=np.abs(x)  
        if is_y_bigger:                            #the numerator is bigger
            
            if y < 0:
                sign = -1
            else:
                sign = 1
            
            ratio = np.round(y/x)                  #find step size for y
            y_step = int(ratio)
            x_step = 1

            iterations = np.ceil(y/ratio)          #find how many times ratio divides into y

            
            
            if iterations < x:                     #if it's less than total x steps, go x dir first

                loc = [loc[0]+1, loc[1]]           #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                x=x-x_step

                if check_where(loc, grid, 'right'):
                    return curve

            else:                                  #otherwise, continue as normal

                x=x
                
            y_count = np.abs(y)
            x_count = np.abs(x)

            while y_count>0 and x_count>0:                 #we'll stop when we've gone the total distance
                for i in range(sign*y_step):        #do y steps up, show each one
                    loc = [loc[0] , loc[1]+sign*1]
                    here = grid[loc[1], loc[0]]

                    pt[loc[1], loc[0]] = 10
                    mask = pt > 0
                    g1 = np.empty(mask.shape)
                    g1.fill(np.nan)
                    g1[mask == True] = 1

                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.savefig('dec_'+str(count)+'.png')
                        count+=1
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.savefig('dec_'+str(count)+'.png')
                        count+=1
                        plt.show()
                    
                    unit = check_type(here)
            
                    curve.append(float(here/unit))
                    time.sleep(1/2)
                    
                    if y > 0:

                        if check_where(loc, grid, 'up'):
                            return curve
                    else:
                        
                        if check_where(loc, grid, 'down'):
                            return curve

                loc = [loc[0]+1, loc[1]]      #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()
                
                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y_count=y_count-y_step                     #update total y left to go
                x_count=x_count-x_step                     #update total x left to go
                if y_count < y_step:                 #if total remaining y is less than a total step, then the step becomes the remainder
                    y_step = y_count

                if check_where(loc, grid, 'right'):
                    return curve


        else:                        #the denominator is bigger
            ratio = np.round(x/y)                  #find step size for x
            y_step = 1
            x_step = int(ratio)

            iterations = np.ceil(x/ratio)          #find how many times ratio divides into x
            
            if x < 0:
                sign = -1
            else:
                sign = 1

            if iterations < y:                     #if it's less than total y steps, go y dir first

                loc = [loc[0], loc[1]+1]           #do 1 step up, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()

                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y=y-y_step

                if check_where(loc, grid, 'up'):
                    return curve

            else:                                  #otherwise, continue as normal

                y=y
                
            y_count = np.abs(y)
            x_count = np.abs(x)

            while y_count>0 and x_count>0:                 #we'll stop when we've gone the total distance
                for i in range(sign*x_step):        #do x steps right, show each one
                    loc = [loc[0]+1 , loc[1]]
                    here = grid[loc[1], loc[0]]

                    pt[loc[1], loc[0]] = 10
                    mask = pt > 0
                    g1 = np.empty(mask.shape)
                    g1.fill(np.nan)
                    g1[mask == True] = 1

                    if w:
                        fig = plt.figure(figsize=(5,5))
                        ax = fig.add_subplot(111, projection=w)
                        im1 = ax.imshow(grid, cmap=c)
                        im2 = ax.imshow(g1*pt, cmap='Greys')

                        ra = ax.coords[0]
                        ra.set_major_formatter('hh:mm:ss')
                        ra.set_separator(':')
                        dec = ax.coords[1]
                        ax.set_xlabel('RA (J2000)')
                        ax.set_ylabel('DEC (J2000)')
                        plt.savefig('dec_'+str(count)+'.png')
                        count+=1
                        plt.show()

                    else:
                        plt.figure(figsize=(5, 5))
                        plt.pcolormesh(grid)
                        plt.pcolormesh(g1*pt, cmap='Greys')
                        plt.savefig('dec_'+str(count)+'.png')
                        count+=1
                        plt.show()
                    
                    unit = check_type(here)
            
                    curve.append(float(here/unit))
                    time.sleep(1/2)

                    if check_where(loc, grid, 'right'):
                        return curve

                loc = [loc[0], loc[1]+sign*1]      #do 1 step right, show it
                here = grid[loc[1], loc[0]]
                pt[loc[1], loc[0]] = 10
                mask = pt > 0
                g1 = np.empty(mask.shape)
                g1.fill(np.nan)
                g1[mask == True] = 1

                if w:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_subplot(111, projection=w)
                    im1 = ax.imshow(grid, cmap=c)
                    im2 = ax.imshow(g1*pt, cmap='Greys')

                    ra = ax.coords[0]
                    ra.set_major_formatter('hh:mm:ss')
                    ra.set_separator(':')
                    dec = ax.coords[1]
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()

                else:
                    plt.figure(figsize=(5, 5))
                    plt.pcolormesh(grid)
                    plt.pcolormesh(g1*pt, cmap='Greys')
                    plt.savefig('dec_'+str(count)+'.png')
                    count+=1
                    plt.show()
                
                unit = check_type(here)
            
                curve.append(float(here/unit))
                time.sleep(1/2)

                y_count=y_count-y_step                     #update total y left to go
                x_count=x_count-x_step                     #update total x left to go
                if x_count < x_step:                 #if total remaining x is less than a total step, then the step becomes the remainder
                    x_step = x_count

                if check_where(loc, grid, 'up'):
                    return curve
   
    
    
    return curve



    
def fill_nan(data):
    """
    For a given data cube, find all columns with nans, and fill with noise based off the rest of the data
    
    Fills for mean, std based off of entire set, but an improvment could be to fill blank space only based off of data where there is no missible information.
    """
    
    mean_test = np.ravel(np.mean(data, 0))
    std_test = np.ravel(np.std(data, 0))
    mean = np.mean(mean_test[TtF(np.isnan(mean_test))])
    std = np.mean(std_test[TtF(np.isnan(std_test))])

    mean=0


    for i in range(data[0].shape[1]):
        for j in range(data[0].shape[0]):
            if np.isnan(data[0][j][i]):
                noise = np.random.normal(mean,std,80)
                data.T[i][j] = noise
            else:
                continue
    return None

def check_type(var):
    """
    Check if data type is an astropy quantity
    """
    
    if type(var) == astropy.units.quantity.Quantity:
        return var.unit
    else:
        return 1
            
            
            
            
            
            
            


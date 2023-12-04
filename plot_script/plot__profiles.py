from math import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib.font_manager as font_manager
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset, zoomed_inset_axes)


# load data
table_Analytical__TM   = '../000001_TM.dat'


# settings
FileOut = 'fig__profiles'

f, ax = plt.subplots( 3, 1, sharex=False, sharey=False )
f.subplots_adjust( hspace=0.05, wspace=0.3 )
f.set_size_inches( 15, 15 )



# pressure
Analytical__TM_pres        = np.loadtxt( table_Analytical__TM,  usecols=(3),  unpack=True )

# density
Analytical__TM_dens        = np.loadtxt( table_Analytical__TM,  usecols=(1),  unpack=True )

# 4-velocity
Analytical__TM_4vel        = np.loadtxt( table_Analytical__TM,  usecols=(2),  unpack=True )

# position
position                   = np.loadtxt( table_Analytical__TM,  usecols=(0),  unpack=True )

Color_TM   = 'blue'


# plot pressure, density, velocity
# ============================================================
ax[0].plot(position , Analytical__TM_pres  ,'-' , color=Color_TM   ,lw=2, label='Analytical (TM)'          )
ax[0].set_ylabel( 'Pressure',       fontsize=20,  fontweight='bold' )
ax[0].tick_params( which='both', tick2On=True, direction='in', labelsize=20, pad=10  )
ax[0].set_xlim(min(position), max(position))
ax[0].get_xaxis().set_ticks([])
ax[0].tick_params( which='both', tick2On=True, direction='in', labelsize=20 , pad=10 )
ax[0].get_xaxis().set_ticks([])

ax[1].plot(position, Analytical__TM_4vel,  '-' , color=Color_TM   ,lw=2                                                        , label='Exact (TM EoS)'      )
ax[1].set_ylabel( '4-velocity',       fontsize=20,  fontweight='bold'             )
ax[1].tick_params( which='both', tick2On=True, direction='in', labelsize=20, pad=10  )
ax[1].set_xlim(min(position), max(position))
ax[1].get_xaxis().set_ticks([])

ax[2].plot(position, Analytical__TM_dens,     '-', color=Color_TM    , lw=2                                           , label='Exact (TM EoS)'       )
ax[2].set_ylabel( r'$\rho$',       fontsize=20,  fontweight='bold' )
ax[2].tick_params( which='both', tick2On=True, direction='in', labelsize=20, top=False, pad=10  )
ax[2].set_xlabel('x', fontsize=20)
ax[2].set_xlim(min(position), max(position))




# save/show figure
# ============================================================
plt.savefig( FileOut+".png", bbox_inches='tight', pad_inches=0.05, format='png' )

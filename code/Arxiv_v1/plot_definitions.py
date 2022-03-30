def make_raw_arrays(filename) :
    df = pd.read_csv(filename)
    
    x = np.array(df['False negative test rate'])
    y = np.array(df['Tracing efficiency'])
    z = np.array(df['Reff '])   
    return x,y,z

def make_cbar_arrays(x,y,z,smooth_horizontally = True) :
    x=np.unique(x)
    y=np.unique(y)


    X,Y = np.meshgrid(list((x)),list((y)))
    Z=z.reshape(len(y),len(x))
    if (smooth_horizontally == True) :
        for i in range (np.shape(Z)[0]):
            Z[i,:] = gaussian_filter(Z[i,:],0)
    return X,Y,Z


def get_contours(x,y,z,delays) :
    
    contours = {}
    for key in delays : 
        contours[key] = {'x':[],'y':[],'z':[]}
        for i in range (len(y)) :
            contours[key]['x'].append(x[i])
            contours[key]['y'].append(y[i])
            contours[key]['z'].append(z[i]-delays[key]['ydict'][y[i]])
    return contours
        

def plot_heatmap(X,Y,Z,delays,interpolate = False,figname=None,heatmap='viridis',R_contours=[1], smooth_horizontally = True) :
    fontsize_medium = 14
    fontsize_small = 12
    
    
    
    plt.figure()
    plt.pcolormesh(X,Y,gaussian_filter(Z,0),cmap=heatmap)
    
    
    cbar = plt.colorbar()

 
    
    for key in delays.keys() :
        ZZ = copy.copy(delays[key]['Z'])
      
        ZZ = gaussian_filter(delays[key]['Z'],1.)   
                
                
        ct=plt.contour(delays[key]['X'],delays[key]['Y'],ZZ,[0],cmap=cm.gray) 
        plt.clabel(ct,ct.levels, inline=True, fmt='%s days'%key)
    plt.xlabel('False negative rate',fontsize=fontsize_medium)
    plt.ylabel('Tracing efficiency',fontsize=fontsize_medium)
    cbar.set_label(r'$R_{\rm eff}$',fontsize=fontsize_small)
    
    
    if (figname != None) :
        plt.savefig('Figures/%s.png'%figname,dpi=400)
    
    plt.show()
    
    
def plot_heatmap_plusminus(X,Y,Z,delays,delay,interpolate = False,figname=None,heatmap='bwr',R_contours=[1], smooth_horizontally = True) :
    fontsize_medium = 14
    fontsize_small = 12
    
    
    
    plt.figure()

 
    
    #for key in delays.keys() :
    ZZ = copy.copy(delays[delay]['Z'])                
    ZZ = gaussian_filter(delays[delay]['Z'],1.)   
                
                
    
    
    plt.pcolormesh(delays[delay]['X'],delays[delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap)
    
    plt.clim(-1, 1)
    cbar = plt.colorbar()

    plt.xlabel('False negative rate',fontsize=fontsize_medium)
    plt.ylabel('Tracing efficiency',fontsize=fontsize_medium)
    cbar.set_label(r'$\Delta R_{\rm eff}$',fontsize=fontsize_small)
    
    plt.text(0.80,.05,"Wait: "+"%i"%delay+" days")
    if (delay==4) :
        plt.text(0.25,0.95,"Fast better",color='w')
        plt.text(0.70,0.95,"Accurate better",color='w')


    
    
    plt.xlim([0.20,1])
    
    if (figname != None) :
        plt.savefig('Figures/%s.png'%figname,dpi=400)
    plt.figure(figsize=(8,8/20))
    plt.pcolormesh(delays[delay]['X'],delays[delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap)
    plt.xlim([0.20,1])
    plt.ylim([0.80,0.80])
    plt.yticks([])
    plt.clim(-.75, .75)
    plt.savefig('Figures/bar_%s.png'%delay)
    
    plt.show()    
    
    
def make_fig2(X,Y,Z,delays,interpolate = False,figname=None,heatmap='bwr',R_contours=[1], smooth_horizontally = True) :


    # MAKE DEFINITIONS FOR PLOT

    # Getting numbering of panels in the correct position.....

    divide_h = 40

    def h_text(xlim,width='narrow') :

        if ( width == 'wide') :
            return (xlim[1] - (xlim[1]-xlim[0])/divide_h)

        elif ( width == 'medium') :
            return (xlim[1] - (xlim[1]-xlim[0])/divide_h*1.5)    

        else :     

            return (xlim[1] - (xlim[1]-xlim[0])/divide_h*3)
    def v_text(ylim,width='narrow') :

        if ( width == 'wide') :
            return(ylim[1] - (ylim[1]-ylim[0])/7)#10)

        else : 
            return(ylim[1] - (ylim[1]-ylim[0])/5)#7)

    # convert cm to inches
    def cm_to_inch(cm) :
        return 0.3937007874*cm

    # Plot this percentile in main frames..
    percentile = 50
    N_samples = 2



    # Import packages for plots
    from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                      mark_inset)
    from matplotlib import cm


    # Fig 2
    figsize2 = (cm_to_inch(17.8),cm_to_inch(17.8/3*2))

    fig = plt.figure(figsize=figsize2)

    # General: 
    plot_dimension = (24,3)
    alpha50 = .6#.25#.50
    alpha90 = alpha50/2#.25
    hm_scattersize = 2#1.#1.5#2#5
    hm_small_scattersize = 1.#1.5#2#5
    peak_cmap =  ['viridis','magma','cividis','plasma','inferno'][0]
    hm_range = [-1.25,1.25]
    hm_fixy = 0.8
    
    outside_color = 'C4'
    outside_marker = 'v'
    outside_alpha=0.8


    datamarkersize=1
    FT_linewidth = 0.5#1
    FT90_linewidth = 0.5
    alphaFT90 = 0.5

    # INSET NUMBERING
    numberingfontsize = 14
    halignment = 'right'
    valignment='center'

    # FONTS
    EXTREMELY_SMALL_SIZE = 5#4.2
    EVEN_SMALLER_SIZE = 6.3#5.5
    SMALL_SIZE = 7
    MEDIUM_SIZE = 8#10
    BIGGER_SIZE = 12
    
    # LABELS AND TEXT
    wait_xpos = 0.675-0.10

    # GENERAL SETTINGS
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=EVEN_SMALLER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # BOXPLOT PLOT SETTINGS
    BP_xticks = (-50,0,50,100,150)#(0,400,800)
    BP_ylim = [-10,850]#[-200,8800]

    hide_yticklabels = True
    hide_xticklabels = False

    show_inset_labels = True
    inset_tickwidth = .4
    inset_ticklength = 2
    inset_pad =1

    horizontal_space_between_subplots = 5.#1.
    vertical_space_between_subplots = .4#.2

    fontsize_medium = 14
    fontsize_small = 12
    
    
    gauss_smooth = .25
    
    
    
    # ---
    # TOP: Bars
    # ---

    # define
    ax0 = plt.subplot2grid(plot_dimension,(4,0),rowspan=2,colspan=1)
    
    plot_delay = 2
    ZZ = copy.copy(delays[plot_delay]['Z'])                
    ZZ = gaussian_filter(delays[plot_delay]['Z'],gauss_smooth)     
    
    
    ax0.pcolormesh(delays[plot_delay]['X'],delays[plot_delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap,vmin=hm_range[0],vmax=hm_range[1])
    ax0.set_ylim([hm_fixy,hm_fixy])
    ax0.set_yticks([])
    ax0.set_yticklabels([])

    
    ax0.set_xlabel('False negative rate',fontsize=MEDIUM_SIZE,labelpad=inset_pad)
    ax0.xaxis.set_tick_params(pad=inset_pad)
    ax0.yaxis.set_tick_params(pad=inset_pad)

    xlim = ax0.get_xlim()
    ylim = ax0.get_ylim()
    
    
    ax0.annotate('Breaking point', xy=(0.09, 1), xytext=(0., 20),xycoords="axes fraction",textcoords="offset points",
                 va="center",ha="center",fontsize=EVEN_SMALLER_SIZE,
            arrowprops=dict(arrowstyle='-|>',facecolor='black'),
            )

    ax0.text(1.00,.85,"Accurate better",fontsize=EXTREMELY_SMALL_SIZE,ha="right")    
        
    ax0.text(-0.10,1.02,r"A",fontweight='bold',fontsize=BIGGER_SIZE,ha="left",va="top")
    
    
    # Bar, delay 3
    ax1 = plt.subplot2grid(plot_dimension,(12-1,0),rowspan=2,colspan=1)
    
    plot_delay = 3
    ZZ = copy.copy(delays[plot_delay]['Z'])                
    ZZ = gaussian_filter(delays[plot_delay]['Z'],gauss_smooth)     
    
    
    ax1.pcolormesh(delays[plot_delay]['X'],delays[plot_delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap,vmin=hm_range[0],vmax=hm_range[1])
    ax1.set_ylim([hm_fixy,hm_fixy])
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    
    ax1.set_xlabel('False negative rate',labelpad=inset_pad,fontsize=MEDIUM_SIZE)
    ax1.xaxis.set_tick_params(pad=inset_pad)
    ax1.yaxis.set_tick_params(pad=inset_pad)

    ax1.annotate('', xy=(0.37, 1), xytext=(0., 20),xycoords="axes fraction",textcoords="offset points",
                 va="center",ha="center",
            arrowprops=dict(arrowstyle='-|>',facecolor='black'),
            )    
    ax1.text(0.0,.85,"Fast better",fontsize=EXTREMELY_SMALL_SIZE,ha="left")    
    ax1.text(1.,.85,"Accurate better",fontsize=EXTREMELY_SMALL_SIZE,ha="right")    
    
    # Bar, delay 4
    ax2 = plt.subplot2grid(plot_dimension,(21-3,0),rowspan=2,colspan=1)
    
    plot_delay = 4
    ZZ = copy.copy(delays[plot_delay]['Z'])                
    ZZ = gaussian_filter(delays[plot_delay]['Z'],gauss_smooth)     
    
    
    ax2.pcolormesh(delays[plot_delay]['X'],delays[plot_delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap,vmin=hm_range[0],vmax=hm_range[1])
    ax2.set_ylim([hm_fixy,hm_fixy])
    ax2.set_yticks([])
    
    ax2.set_xlabel('False negative rate',fontsize=MEDIUM_SIZE,labelpad=inset_pad)
    ax2.xaxis.set_tick_params(pad=inset_pad)
    ax2.yaxis.set_tick_params(pad=inset_pad)

    
    ax2.annotate('', xy=(0.54, 1), xytext=(0., 20),xycoords="axes fraction",textcoords="offset points",
                 va="center",ha="center",
            arrowprops=dict(arrowstyle='-|>',facecolor='black'),
            )        
    ax2.text(0.0,.85,"Fast better",fontsize=EXTREMELY_SMALL_SIZE,ha="left")    
    ax2.text(1.,.85,"Accurate better",fontsize=EXTREMELY_SMALL_SIZE,ha="right")    
        
    '''
    MIDDLE 1 : 
    '''

    # ----
    # AX01: Delay 2
    # ----
    plot_delay = 2
    ax01 = plt.subplot2grid(plot_dimension,(0,1),rowspan=11,colspan=1)

    #for key in delays.keys() :
    ZZ = copy.copy(delays[plot_delay]['Z']) 
    ZZ = gaussian_filter(delays[plot_delay]['Z'],gauss_smooth)       
    
    ax01.pcolormesh(delays[plot_delay]['X'],delays[plot_delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap,vmin=hm_range[0],vmax=hm_range[1])
    ax01.plot(np.arange(0,1,0.02),0.8*np.ones(50),'k--',linewidth=0.5,alpha=0.25)

    plt.ylabel('Tracing efficiency',fontsize=MEDIUM_SIZE,labelpad=inset_pad)
    plt.text(wait_xpos,.05,"Wait: "+"%i"%plot_delay+" days",fontsize=SMALL_SIZE)
       

    ax01.set_xticklabels([])
    ax01.text(0.10-.35,1.0,r"B",fontweight='bold',fontsize=BIGGER_SIZE,ha="left",va="top")
    
    
    # Labels and other settings


    # ----
    # AX02: Delay 3
    # ----

    plot_delay = 3
    ax02 = plt.subplot2grid(plot_dimension,(0,2),rowspan=11,colspan=1)

    ZZ = copy.copy(delays[plot_delay]['Z']) 
    ZZ = gaussian_filter(delays[plot_delay]['Z'],gauss_smooth)       
    
    im1 = ax02.pcolormesh(delays[plot_delay]['X'],delays[plot_delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap,vmin=hm_range[0],vmax=hm_range[1])
    ax02.plot(np.arange(0,1,0.02),0.8*np.ones(50),'k--',linewidth=0.5,alpha=0.25)

    ax02.text(wait_xpos,.05,"Wait: "+"%i"%plot_delay+" days",fontsize=SMALL_SIZE)
    ax02.text(0.10-0.25,1.0,r"C",fontweight='bold',fontsize=BIGGER_SIZE,ha="left",va="top")

    
    ax02.set_xticklabels([])
    ax02.set_yticklabels([])
    
    
    offset = 0.438
    cbar_ax = fig.add_axes([.91, offset+0.1, 0.02, 0.343])
    
    fig.colorbar(im1, cax=cbar_ax,label=r'$\Delta R_{\rm eff}$')        
    


    '''
    MIDDLE 2 :
    '''


    # ----
    # AX11: Delay 4
    # ----

    
    plot_delay = 4
    ax11 = plt.subplot2grid(plot_dimension,(13,1),rowspan=11,colspan=1)

    ZZ = copy.copy(delays[plot_delay]['Z']) 
    ZZ = gaussian_filter(delays[plot_delay]['Z'],gauss_smooth)       
    
    ax11.pcolormesh(delays[plot_delay]['X'],delays[plot_delay]['Y'],gaussian_filter(ZZ,0),cmap=heatmap,vmin=hm_range[0],vmax=hm_range[1])
    ax11.plot(np.arange(0,1,0.02),0.8*np.ones(50),'k--',linewidth=0.5,alpha=0.25)

    ax11.set_xlabel('False negative rate',fontsize=MEDIUM_SIZE,labelpad = inset_pad)
    ax11.set_ylabel('Tracing efficiency',fontsize=MEDIUM_SIZE,labelpad=inset_pad)
    ax11.text(wait_xpos,.05,"Wait: "+"%i"%plot_delay+" days",fontsize=SMALL_SIZE)

    ax11.text(0.10-0.35,1.0,r"D",fontweight='bold',fontsize=BIGGER_SIZE,ha="left",va="top")
    
   
    # Labels and other settings

    # ----
    # AX12 : Colormap Reff
    # ----
    ax12 = plt.subplot2grid(plot_dimension,(13,2),rowspan=11,colspan=1)

    im = plt.pcolormesh(X,Y,gaussian_filter(Z,0),cmap='viridis')
    

 
    
    for key in delays.keys() :
        ZZ = copy.copy(delays[key]['Z'])
        ZZ = gaussian_filter(delays[key]['Z'],gauss_smooth)   
  
        ct=plt.contour(delays[key]['X'],delays[key]['Y'],ZZ,[0],cmap=cm.gray) 
        plt.clabel(ct,ct.levels, inline=True, fmt='%s days'%key)

    ax12.set_xlabel('False negative rate',fontsize=MEDIUM_SIZE)   
    
    # Labels and other settings

    ax12.set_yticklabels([])
    


    cbar_ax = fig.add_axes([.91, 0.125, 0.02, 0.343])
    fig.colorbar(im, cax=cbar_ax,label=r'$R_{\rm eff}^{\rm rapid}$')    

    ax12.text(0.10-.25,1.0,r"E",fontweight='bold',fontsize=BIGGER_SIZE,ha="left",va="top")
 

    box = ax0.get_position()
    box.x0 = box.x0  -0.05
    box.x1 = box.x1 - 0.05
    ax0.set_position(box)    
    
    box = ax1.get_position()
    box.x0 = box.x0  -0.05
    box.x1 = box.x1 - 0.05
    ax1.set_position(box)  
    
    box = ax2.get_position()
    box.x0 = box.x0  -0.05
    box.x1 = box.x1 - 0.05
    ax2.set_position(box)  
    
    
    plt.savefig('Figures/Fig2.png',dpi=400)
    plt.savefig('Figures/Fig2.pdf',dpi=400)
    plt.savefig('Figures/Fig2.svg',dpi=400)
    

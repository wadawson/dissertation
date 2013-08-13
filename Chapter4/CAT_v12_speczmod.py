'''
Cluster Analysis Tool

This program calculates the redshift and velocity dispersion for a galaxy
cluster, based on an input spectrocopic catalog containing ra, dec, & redshift
for a population of cluster galaxies. The biweight statistic and bias corrected
percent confidence limit (Beers et al. 1990) are used.

As of version 1.0 it includes the number density mapping function numberdensity
which creates a 2D histogram based map of the galaxies and creates a fits file
of that histogram complete with wcs

Currently the three main functions are:
zhist - Creates a redshift histogram for all the galaxies with redshifts within
        'zrange' and specified physical 'radii'.
        
zVdisp - Estimates the redshift and velocity dispersion of the specified cluster
         using the galaxies within 'vwidth' of the cluster redshift and within
         the specified physical 'radii'.
         
         
version 0.0
-- contains zhist, and zVdisp functions
version 1.0
-- added numberdensity
version 1.1
-- modified zVdisp such that it now outputs pickled arrays of the Location and
   scale bootstrap analysis
version 1.2
-- added a bootstrap option to the numberdensity map function
version speczmod
-- modified to incorporate spec-z + photo-z info.  This will remove photo-z info
   of any galaxy with spec-z infor, and will fix all spec-z galaxies during the
   bootstrap resampling.
'''
from __future__ import division
import numpy
import pylab
import tools
import sys
import cosmo


def zhist(catalog, colnames, cluster_coord, zbinwidth, zrange, radii, filename,
          h = 0.7, Om = 0.3,Ol = 0.7):
    '''
    Creates a redshift histogram with bins of width 'zbinwidth', for all the 
    galaxies with redshifts within 'zrange' and specified physical 'radii'.
    
    Input:
    catalog = [string] a catalog with ra, dec, and z columns. It must have a
        ttype header either 0 or 1 indexed, e.g.:
        #ttype1 = 'ra'
        #ttype2 = 'dec'
        #ttype3 = 'z'
    colnames = [(string,string,string)] the ttype names of (ra, dec, z)
    cluster_coord = [(float,float); units:(deg,deg)] the cluster center (ra,dec)
    zbinwidth = [float; units:km/s] width of the redshift and velocity 
        dispersion histogram bins
    zrange = [(float,float)] redshift range to consider (z_min, z_max)
    radii = [List of floats; units:Mpc] List of the radial bins about cluster
        center to analyze. A separate z and velocity dispersion will be 
        estimated using the cluster redshift galaxies within each radii.
    filename = [string] name of the output histogram png
    h = [float] hubble scale factor, i.e. H = 100km/s/h
    Om = [float] matter energy density
    Ol = [float] dark energy density
    
    Output:

    '''
    c = 3.e5 #speed of light in km/s

    ### CALCULATIONS ################
    # Read in the catalog data and header
    cat = tools.readcatalog(catalog)
    key = tools.readheader(catalog)
    racol = key[colnames[0]]
    deccol = key[colnames[1]]
    zcol = key[colnames[2]]
    # filter the catalog for the user specified redshift range
    cat = tools.filtercat(cat,zcol,zrange[0],zrange[1])
    ## Create redshift histograms ########
    # create the histogram reshift bin array
    zmin = numpy.min(cat[:,zcol])
    zmax = numpy.max(cat[:,zcol])
    # consirvatively estimate the number of bins
    Nest = c*(zmax-zmin)/(1+zmin)/zbinwidth
    zbin = numpy.zeros(Nest+50)
    zbin[0] = zmin
    for i in numpy.arange(Nest):
        zbin[i+1]= (zbinwidth/c*(1+zbin[i]/2)+zbin[i])/(1-zbinwidth/c/2)
        if zbin[i+1] > zmax:
            zbin = zbin[:i+2]
            break
    #create the histogram plot for each of the radii
    #loop through the catalog and measure the angular separation of each
    #object from the center of the cluster
    ang = numpy.zeros(numpy.shape(cat[:,0])[0])
    for i in numpy.arange(numpy.shape(cat[:,0])[0]):
        ang[i] = tools.angdist(cluster_coord[0],cluster_coord[1],cat[i,racol],cat[i,deccol])*60
    # since the angular diamenter distance will change from redshift bin to 
    # redshift bin it is necessary to create a unique set of radial filters for 
    # each redshift bin
    counts = numpy.zeros((numpy.size(zbin)-1,numpy.size(radii)))
    for i in numpy.arange(numpy.size(zbin)-1):
        zavg = (zbin[i]+zbin[i+1])/2
        for j in numpy.arange(numpy.size(radii)):
            maxang = radii[j]/cosmo.ProjectedLength(zavg,h=h,Om=Om,Ol=Ol)
            mask_ang = ang <= maxang
            cat_flt = cat[mask_ang,:]
            mask_z = numpy.logical_and(cat_flt[:,zcol]>=zbin[i],
                                       cat_flt[:,zcol]<zbin[i+1])
            counts[i,j] = numpy.sum(mask_z)
    fig1 = pylab.figure()
    for i in numpy.arange(numpy.size(radii)):
        pylab.step(zbin[:-1],counts[:,i],label='R<{0}Mpc'.format(radii[i]),
                   where='post')
    pylab.legend(loc=0)
    pylab.title('Cluster at Ra={0:0.3f} dec={1:0.2f}'.format(cluster_coord[0],
                                                             cluster_coord[1]))
    pylab.xlabel('Redshift Bins of {0}(1+z)km/s'.format(zbinwidth))
    pylab.ylabel('$N_{galaxy}$')
    pylab.savefig(filename)
    pylab.show()

def zVdisp(catalog, colnames, cluster_coord, zest, vwidth, radii, prefix,
           Nboot = 10000,N_sigma = 1, zbinwidth = 250, h = 0.7, Om = 0.3,
           Ol = 0.7):
    '''
    Estimates the redshift and velocity dispersion of the specified cluster
    using the galaxies within 'vwidth' of the cluster redshift and within the
    specified physical 'radii'.
    
    Input:
    catalog = [string] a catalog with ra, dec, and z columns. It must have a
        ttype header either 0 or 1 indexed, e.g.:
        #ttype1 = 'ra'
        #ttype2 = 'dec'
        #ttype3 = 'z'
    colnames = [(string,string,string)] the ttype names of (ra, dec, z)
    cluster_coord = [(float,float); units:(deg,deg)] the cluster center (ra,dec)
    zest = [float] initial guess at the cluster redshift
    vwidth = [float; units:km/s] galaxies within +/-'vwidth' of the cluster
        redshift will be used in the estimate of the cluster redshift and
        velocity dispersion. This should be at least three times the expected
        velocity dispersion, so ~3000. Note that the biweight will down weight
        the tails.
    radii = [List of floats; units:Mpc] List of the radial bins about cluster
        center to analyze. A separate z and velocity dispersion will be 
        estimated using the cluster redshift galaxies within each radii.
    prefix = [string] prefix to append to each of the output filenames
    Nboot = [int] number of bootstrap samples to generate when estimating the
        confidence intervals for z and velocity dispersion.
    N_sigma = [int] the number of sigma to report for the confidence limits e.g.
        for 95% confidence limits N_sigma=2
    zbinwidth = [float; units:km/s] width of the redshift and velocity 
        dispersion histogram bins
    h = [float] hubble scale factor, i.e. H = 100km/s/h
    Om = [float] matter energy density
    Ol = [float] dark energy density
    
    Output:
    
    '''
    import pickle
    from scipy.stats import norm
    from scipy.special import erf
    
    c = 3.e5 #speed of light in km/s
    
    ### Functions ################
    
    def biweightLoc(z,c=6):
        M = numpy.median(z)
        MAD = numpy.median(numpy.abs(z))
        u = (z-M)/(c*MAD)
        mask_u = numpy.abs(u) < 1
        z = z[mask_u]
        u = u[mask_u]
        Cbi = M + numpy.inner(z-M,(1-u**2)**2)/numpy.sum((1-u**2)**2)
        return Cbi
    
    def biweightScale(z,c=9):
        n = numpy.size(z)
        M = numpy.median(z)
        MAD = numpy.median(numpy.abs(z))
        u = (z-M)/(c*MAD)
        mask_u = numpy.abs(u) < 1
        z = z[mask_u]
        u = u[mask_u]
        Sbi = n**(0.5)*numpy.inner((z-M)**2,(1-u**2)**4)**(0.5)/numpy.abs(numpy.inner(1-u**2,1-5*u**2))
        return Sbi
    
    def bcpcl(T,T_p,N_sigma):
        '''
        Calculates the bias corrected percent confidence limits.
        -- Suppose that we have observed data (y1, y2, ..., yn) and use it to estimate a population parameter Q (e.g. Q could be the true mean of the entire population).
        -- T is a statistic that estimates Q. For example T could be an estimate of the true mean by calculating the mean of  (y1, y2, ..., yn).
        -- Suppose that we create m bootstrap samples (y_p_1j, y_p_2j, ...,j_p_nj) from observed sample  (y1, y2, ..., yn), where j is the jth bootstrap sample.
        -- Then T_p_j is the jth bootstrap observation of T.  For example this could be the mean of (y_p_1j, y_p_2j, ...,j_p_nj).
        
        T = [float] e.g. biweight Location for (y1, y2, ..., yn)
        T_p = [vector array] biwieght Locations for the bootstrap samples
        N_sigma = the number of sigma to report the confidence limits for
            e.g. for 95% confidence limits N_sigma=2
        Return (lower, upper) confidence limits
        '''
        #Percentile confidence interval is defined as 100%(1-a), thus for 1sigma a=0.32
        a = 1-erf(N_sigma/numpy.sqrt(2))
        #order the bootstrap sample values smallest to largest
        index = numpy.argsort(T_p)
        T_p = T_p[index]
        #Number of bootstrap samples
        m = numpy.size(T_p)        
        #Calculate the bias correction term
        mask = T_p < T
        z_0 = norm.ppf(numpy.sum(mask)/m)
        #Calculate the a1 and a2 values
        a1 = norm.cdf(2*z_0+norm.ppf(a/2))
        a2 = norm.cdf(2*z_0+norm.ppf(1-a/2))
        #Calculate the lower and upper indicies of lower and upper confidence intervals
        id_L = numpy.int(m*a1)-1
        id_U = numpy.int(m*a2)
        #Find the lower an upper confidence values
        T_L = T_p[id_L]
        T_U = T_p[id_U]
        return (T_L, T_U)
    
    def trimcat(cat,z,zcol,vwidth,radius,angcol):
        '''
        cat = [array]
        z = [float] cluster redshift
        zcol = [int] the redshift column number (zero indexed)
        vwidth = [float; units:km/s] the generous +/- velocity bounds of the 
                 cluster
        radius = [float; units:Mpc] only consider galaxies closer that radius
        radiuscol = [int] the angular separation (arcmin) column number (zero
                    indexed)
        '''
        #Perform the redshift cut
        zmax = (z+vwidth/c*(1+z/2))/(1-vwidth/(2*c))
        zmin = (z-vwidth/c*(1+z/2))/(1+vwidth/(2*c))
        cat =  tools.filtercat(cat,zcol,zmin,zmax,verbose=False)
        #Perfrom the physical separation cut
        angmax = radius/cosmo.ProjectedLength(z)
        return tools.filtercat(cat,angcol,None,angmax,verbose=False)
    
    ### Calculations #################3
    
    # Read in the catalog data and header
    cat = tools.readcatalog(catalog)
    key = tools.readheader(catalog)
    racol = key[colnames[0]]
    deccol = key[colnames[1]]
    zcol = key[colnames[2]]
    
    # Create the blank output arrays
    Nrad = numpy.size(radii)
    Cbi_array = numpy.zeros(Nrad)
    Sbi_array = numpy.zeros(Nrad)
    Cbi_p_array = numpy.zeros((Nrad,Nboot))
    Sbi_p_array = numpy.zeros((Nrad,Nboot))
    Cbi_L_array = numpy.zeros(Nrad)
    Cbi_U_array = numpy.zeros(Nrad)
    Sbi_L_array = numpy.zeros(Nrad)
    Sbi_U_array = numpy.zeros(Nrad)
    
    # Calculate angular separation of the galaxies from the cluster in arcmin 
    ang = numpy.zeros(numpy.shape(cat[:,0])[0])
    for i in numpy.arange(numpy.shape(cat[:,0])[0]):
        ang[i] = tools.angdist(cluster_coord[0],cluster_coord[1],cat[i,racol],
                               cat[i,deccol])*60
    # Add this data to the catalog
    ang = numpy.reshape(ang,(numpy.size(ang),1))
    cat = numpy.concatenate((cat,ang),axis=1)

    for i in numpy.arange(numpy.size(radii)):
        # Iteratively estimate the cluster redshift
        Cbi = zest
        converg = 99
        j = 0
        while converg > 0.0000001:
            fltcat = trimcat(cat,Cbi,zcol,vwidth,radii[i],-1)
            z = fltcat[:,zcol]
            Cbi_old = Cbi
            Cbi = Cbi_array[i] = biweightLoc(z,c=6)
            converg = (Cbi-Cbi_old)/(Cbi+Cbi_old)*2
            j+=1
            if j == 50:
                print "Cluster central redshift hasn't converged after 50 iterations, exiting."
                sys.exit()
        print 'Cluster central redshift estimate converged to {0:0.5f} after {1} iterations.'.format(Cbi,j)
        
        #Calculate the velocities of the galaxies relative to the cluster center
        v = c*(z-Cbi)/(1+(z+Cbi)/2)
        
        # Estimate the cluster velocity scale (i.e. velocity dispersion)
        Sbi = Sbi_array[i] = biweightScale(v,c=9)
                
        # Perform the bootstrap Bias Corrected Percentile Confidence Intervals
        Cbi_p = numpy.zeros(Nboot)
        Sbi_p = numpy.zeros(Nboot)
        
        for j in numpy.arange(Nboot):
#            if j%100 == 0:
#                print 'Processing bootstrap sample {0} of {1}'.format(j,Nboot)
            # create a random bootstrap sample index array
            n = numpy.size(z)
            b = numpy.random.randint(0,high=n,size=n)
            Cbi_p[j] = biweightLoc(z[b])
            Sbi_p[j] = biweightScale(v[b])
        Cbi_p_array[i,:] = Cbi_p
        Sbi_p_array[i,:] = Sbi_p
        ##piclke the bootstrap samples
        #F = open(prefix+'Cbi_p','w')
        #pickle.dump(Cbi_p,F)
        #F.close()
        #F = open(prefix+'Sbi_p','w')
        #pickle.dump(Sbi_p,F)
        #F.close()
                
        # Calculate the confidence limits
        Cbi_L_array[i],Cbi_U_array[i] = bcpcl(Cbi,Cbi_p,N_sigma)
        Sbi_L_array[i],Sbi_U_array[i] = bcpcl(Sbi,Sbi_p,N_sigma)
        
        # Plot a histogram of the redshifts
        fig = pylab.figure()
        pylab.hist(z,bins=20)
        pylab.title('{0} Cluster Galaxies within {1}Mpc of Ra={2:0.3f} dec={3:0.2f}'.format(numpy.size(z),radii[i],cluster_coord[0],cluster_coord[1]))
        pylab.xlabel('Redshift')
        pylab.ylabel('$N_{galaxy}$')
        name = prefix+'_zhist_Rlt{0}Mpc.png'.format(radii[i])
        pylab.savefig(name)

        
    #Print the results
    print "Cluster redshift results; upper and lower {0}sigma confidence limits are quoted:".format(N_sigma)
    print '\t\tCenter\tLower\tUpper\tUnits'
    for i in numpy.arange(numpy.size(radii)):
        print 'For galaxies with r < {0} Mpc'.format(radii[i])
        print 'Location\t{0:0.5f}\t{1:0.5f}\t{2:0.5f}'.format(Cbi_array[i],Cbi_L_array[i],Cbi_U_array[i])
        print 'Scale\t\t{0:0.0f}\t{1:0.0f}\t{2:0.0f}\tkm/s'.format(Sbi_array[i],Sbi_L_array[i],Sbi_U_array[i])
    
    #pickle the bootstrap samples
    filename = prefix+'_bootstrapLoction.pickle'
    F = open(filename,'w')
    pickle.dump(numpy.transpose(Cbi_p_array),F)
    F.close()

    filename = prefix+'_bootstrapScale.pickle'
    F = open(filename,'w')
    pickle.dump(numpy.transpose(Sbi_p_array),F)
    F.close()

        
    #Plot the bootstrap samples
    fig1 = pylab.figure()
    #create the label names
    labname = []
    for i in numpy.arange(numpy.size(radii)):
        labname.append('R<{0}Mpc'.format(radii[i]))
    pylab.hist(numpy.transpose(Cbi_p_array),bins=20,label=labname)
    pylab.title('Histogram of Cluster Location for the Bootstrap Sample')
    pylab.xlabel('Cluster Biweight Location (i.e. redshift)')
    pylab.ylabel('$N_{bootstrap sample}$')
    pylab.legend()
    filename = prefix+'_bootLocHist'
    pylab.savefig(filename)
    
    fig2 = pylab.figure()
    #create the label names
    labname = []
    for i in numpy.arange(numpy.size(radii)):
        labname.append('R<{0}Mpc'.format(radii[i]))
    pylab.hist(numpy.transpose(Sbi_p_array),bins=20,label=labname)
    pylab.title('Histogram of Cluster Scale for the Bootstrap Sample')
    pylab.xlabel('Cluster Biweight Scale (i.e. velocity dispersion km/s)')
    pylab.ylabel('$N_{bootstrap sample}$')
    pylab.legend()
    filename = prefix+'_bootVdispHist'
    pylab.savefig(filename)
    
    pylab.show()
    
def numberdensity(catalog, colnames, rabin, prefix, rarange=None, decrange=None,
                  zrange=None,N_boot=None,speczname=None,speczclustrange=None,
                  magcol=None,magrange=None):
    '''
    Creates a 2D fits map which is essentially a 2D histogram of the galaxies.
    
    Input:
    catalog = [string] a catalog with ra, dec, and optional z columns. It must
        have a ttype header either 0 or 1 indexed, e.g.:
        #ttype1 = 'ra'
        #ttype2 = 'dec'
        #ttype3 = 'z'
    colnames = [(string,string,string) or (string,string)] the ttype names of
        (ra, dec, z) or (ra, dec), where here z refers to photo-z
    rabin = [integer] the number of bins along the ra axis for the 2D
        histogram, also the number of pixels along the ra axis in the
        fits file. The decbin is calculated such that the 2D bins and pixels
        will be aproximately square.
    prefix = [string] prefix on the names of the output histogram and fits files
    rarange = [(float,float)] {units: (degrees,degrees)} RA range to consider
        (ra_min, ra_max)
    decrange = [(float,float)] {units: (degrees,degrees)} Dec range to consider
        (ra_min, ra_max)
    zrange = [(float,float)] redshift range to consider (z_min, z_max), note
        that colnames must have (ra,dec,z) all specified, where here z refers to
        photo-z
    N_boot = number of bootstrap samples to generate
    speczname = [string] ttype name of the spec-z column
    speczclustrange = [(float,float)] redshift range to consider (zmin, zmax)
        when consider an object as part of the cluster or not
    magcol = [string] ttype name of the magnitude column
    magrange = [(float,float)] min and max magnitudes to consider
        
    Output:
    prefix+'_numberdensity.fits' = a galaxy number density fits file of size xpix by ypix
        If N_boot = None then this fits file just contains one frame, the number density map.  If N_boot != None then a multiframe fits image will be created where the frames are as follows:
        1) number density map
        2) S/N map
        3) standard deviation map
        4 thru N_boot+4) bootstrap sample maps
        
    '''
    import pyfits
    ### CALCULATIONS ################
    # Read in the catalog data and header
    cat = tools.readcatalog(catalog)
    key = tools.readheader(catalog)
    racol = key[colnames[0]]
    deccol = key[colnames[1]]
    magcol = key[magcol]
    if numpy.size(colnames) == 3:
        zcol = key[colnames[2]]
    # Apply the coordinate filters
    if rarange != None:
        ra_min = rarange[0]
        ra_max = rarange[1]
        cat = tools.filtercat(cat,racol,ra_min,ra_max)
    if decrange != None:
        dec_min = decrange[0]
        dec_max = decrange[1]
        cat = tools.filtercat(cat,deccol,dec_min,dec_max)
    if magrange != None:
        cat = tools.filtercat(cat,magcol,magrange[0],magrange[1])
  
    # Extract all the galaxies with spec-z's
    zspeccol = key[speczname]
    mask_spec = cat[:,zspeccol] != -99 #all columns that have spec-z
    cat_specz = cat[mask_spec,:]
    mask_photo = mask_spec == False
    cat_photoz = cat[mask_photo,:]
    # filter on spec-z's
    cat_specz = tools.filtercat(cat_specz,zspeccol,speczclustrange[0],speczclustrange[1])
    # filter on photo-z's
    if zrange != None and numpy.size(colnames) == 3:
        # filter the catalog for the user specified redshift range
        cat_photoz = tools.filtercat(cat_photoz,zcol,zrange[0],zrange[1])
    elif zrange != None and numpy.size(colnames) != 3:
        print 'numberdensity: Warning zrange specified but no z column name specified in colnames input. Exiting'
        sys.exit()
    
    # Calculate the number of bins along the dec axis
    if rarange == None:
        ra_min = numpy.min(cat[:,racol])
        ra_max = numpy.max(cat[:,racol])
    if decrange == None:
        dec_min = numpy.min(cat[:,deccol])
        dec_max = numpy.max(cat[:,deccol])
    dec_mean = (dec_max + dec_min) / 2 * numpy.pi / 180.0 #radians
    decbin = rabin*(dec_max-dec_min)//((ra_max-ra_min)*numpy.cos(dec_mean))
    # create the pixel/bin edge arrays
    ra_binwidth = (ra_max-ra_min)/rabin
    dec_binwidth = (dec_max-dec_min)/decbin
    ra_edge = numpy.arange(ra_min,ra_max+ra_binwidth,ra_binwidth)
    dec_edge = numpy.arange(dec_min,dec_max+dec_binwidth,dec_binwidth)
    
    # Create the blank map_array
    if N_boot != None:
        h = numpy.zeros((3+N_boot,decbin,rabin))
        N = numpy.shape(cat_photoz)[0] #number of rows in filtered catalog
        # Random with replacement bootstrap index array 
        b = numpy.random.randint(0,high=N,size=(N,N_boot))
        print 'numberdensity: will perform bootstrap analysis with {0} random iterations'.format(N_boot)
            
    #create the 2D histogram    
    if N_boot == None:
        cat_flt = numpy.concatenate((cat_specz,cat_photoz),axis=0)
        h, tmp_edge, tmp_edge = numpy.histogram2d(cat_flt[:,deccol],cat_flt[:,racol],bins=(dec_edge,ra_edge))
    else:
        cat_flt = numpy.concatenate((cat_specz,cat_photoz),axis=0)        
        h[0,:,:], tmp_edge, tmp_edge = numpy.histogram2d(cat_flt[:,deccol],cat_flt[:,racol],bins=(dec_edge,ra_edge))
        for i in numpy.arange(N_boot):
            cat_flt = numpy.concatenate((cat_specz,cat_photoz[b[:,i],:]),axis=0)           
            h[3+i,:,:], tmp_edge, tmp_edge = numpy.histogram2d(cat_flt[:,deccol],cat_flt[:,racol],bins=(dec_edge,ra_edge))
        #Create signal/noise and standard deviation maps
        h[2,:,:] = numpy.std(h[3:,:,:],ddof=1)
        h[1,:,:] = h[0,:,:]/h[2,:,:]
    #Create the fits file
    H = h
    hdu = pyfits.PrimaryHDU(numpy.float32(H))
    
    # Calculate the wcs CR**** for the fits header    
    xscale = (ra_max - ra_min) * numpy.cos(dec_mean) / rabin
    yscale = (dec_max - dec_min) / decbin
    crval1 = (ra_max - ra_min) / 2 + ra_min
    crval2 = (dec_max - dec_min) / 2 + dec_min
    crpix1 = rabin / 2
    crpix2 = decbin / 2
    
    # Apply the wcs to the fits header
    hdu.header.update('ctype1', 'RA---TAN')
    hdu.header.update('ctype2', 'DEC--TAN')
    hdu.header.update('crval1', crval1)
    hdu.header.update('crval2', crval2)
    hdu.header.update('crpix1', crpix1)
    hdu.header.update('crpix2', crpix2)
    hdu.header.update('cd1_1',xscale)
    hdu.header.update('cd1_2',0)
    hdu.header.update('cd2_1',0)
    hdu.header.update('cd2_2',yscale)
    filename = prefix+'_numberdensity'
    hdu.writeto(filename+'.fits',clobber=True)

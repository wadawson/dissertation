'''
Input file for the creation of galaxy number density maps for the DLS and 
Med-band photo-z catalogs.
'''
from CAT_v12_speczmod2 import numberdensity

catalog= '/Users/dawson/SkyDrive/Research/Clusters/DLSCL09162953/GalaxyDensity/specz_photoz_based/with2012bRedshifts/MedPerrysClusterCat_2013specz.txt'
ra_ttype = 'alpha'
dec_ttype = 'delta'
spec_z_ttype = 'z_spec'
photo_z_ttype = 'z_mid'
photoz_errorfactor = 0.07
spec_z_sigmarange = 3
photo_z_range = (0.43,0.63)
photo_z_penalty=0.63

rabin=200
N_boot = 10000
magcol = 'R'
magrange = (None,24)
rarange=(138.99750,139.09361)
decrange=(29.79978,29.88314)

# South Subcluster
# note that I have to customize this for the north and south subclusters
z_cluster = 0.53414
vdisp_cluster = 770 #km/s
prefix= 'BVRzMed_24maglim_pzpen_south_10000'
numberdensity(catalog, ra_ttype, dec_ttype, spec_z_ttype, photo_z_ttype, 
              photoz_errorfactor, rabin, prefix, z_cluster, vdisp_cluster, 
              rarange, decrange, N_boot, magcol ,magrange,spec_z_sigmarange,photo_z_range,photo_z_penalty)

# North Subcluster
# note that I have to customize this for the north and south subclusters
z_cluster = 0.53074
vdisp_cluster = 740 #km/s
prefix= 'BVRzMed_24maglim_pzpen_north_10000'
numberdensity(catalog, ra_ttype, dec_ttype, spec_z_ttype, photo_z_ttype, 
              photoz_errorfactor, rabin, prefix, z_cluster, vdisp_cluster, 
              rarange, decrange, N_boot, magcol ,magrange,spec_z_sigmarange,photo_z_range,photo_z_penalty)
#!/usr/local/anaconda/bin/python
"""
 Generate SELFE nudging file for CRITFC CMOP forecast from LiveOcean data
 
 Charles Seaton (cseaton@crtifc.org), Oct 30, 2020
"""
import netCDF4 as nc
import numpy as np
import os, sys
from dateutil import parser
from scipy import interp
from scipy.spatial import cKDTree
from scipy.io import FortranFile
from scipy import interpolate
import argparse

#
#  SELFE file reading functions #
# 
"""
Tuomas Karna 2013-02-15
"""
class vc:
  """
  A class for computingComputes SELFE vertical coordinates.
  """
  def __init__(self, nvrt, nz, h_s, h_c, theta_b=None, theta_f=None, ztot=None,
               sigma=None, cs=None, h0=0.01 ):
    self.nvrt    = nvrt
    self.nz      = nz
    self.h0      = h0
    self.h_s     = h_s
    self.h_c     = h_c
    self.theta_b = theta_b
    self.theta_f = theta_f
    self.ztot = ztot
    self.sigma = sigma
    self.cs = cs
    if cs is None and sigma is not None:
      self.cs = self.compute_cs( sigma )

  @classmethod
  def fromVGridFile( cls, filename, h0=0.01 ) :
  #def fromVGridFile( self, filename, h0=0.01 ) :
    """Parses vgrid file and returns a verticalCoordinates object"""
    print('Reading '+filename+' ...')
    f = open(filename,'r')
    words = f.readline().split()
    nvrt = int( words[0] )
    nz = int( words[1] )
    h_s = float( words[2] )
    print('  nvrt=',nvrt,'nz=',nz,'h_s=',h_s)
    f.readline() # 'Z levels'
    ztot = np.zeros(nvrt)
    sigma = np.zeros(nvrt)
    for i in range(nz) :
      words = f.readline().split()
      #print i,words
      ztot[i] = float(words[1])
    f.readline() # 'S levels'
    words = f.readline().split()
    h_c = float( words[0] )
    theta_b = float( words[1] )
    theta_f = float( words[2] )
    print('  h_c=',h_c, 'theta_b=',theta_b, 'theta_f=',theta_f)
    for i in range(nvrt-nz) :
      words = f.readline().split()
      #print i,words
      sigma[i] = float(words[1])
    return cls(nvrt,nz,h_s,h_c,theta_b,theta_f,ztot,sigma,h0=h0)
    # self.nvrt = nvrt
    # self.nz = nz
    # self.h_s = h_s
    # self.h_c = h_c
    # self.theta_b = theta_b
    # self.tehta_f = theta_f
    # self.ztot = ztot
    # self.sigma = sigma
    # self.h0 = h0
    # return self
    #
    # #return self(nvrt,nz,h_s,h_c,theta_b=theta_b,theta_f=theta_f,ztot=ztot,sigma=sigma,h0=h0)

  def compute_cs(self, sigma) :
    return ( (1-self.theta_b)*np.sinh(self.theta_f*sigma)/np.sinh(self.theta_f) +
             self.theta_b*(  np.tanh(self.theta_f*(sigma+0.5))- np.tanh(self.theta_f*0.5)  ) /2 /np.tanh(self.theta_f*0.5) )

  def computeVerticalCoordinates(self, eta, dp, ztmp=None, kbp2=None) :
    """
    Parameters
    ----------
    eta : array_like (nPoints,)
          water elevation at certain points
    dp  : array_like (nPoints,)
          bathymetry at same points

    Returns
    -------
    Z : array_like (nVert,nPoints)
        z coordinates
    kbp2 : array_like (nPoint,)
        bottom level indices
    iwet : array_like (nPoint,)
        wet mask (per node, isolated wet nodes are not removed like in SELFE)
    """

    #Compute z coordinates (only for 3d)
    #sys.stdout.write('Generating vertical coordinates...')
    #sys.stdout.flush()
    nNodes = len(eta)
    if ztmp is None:
        ztmp = np.zeros( (self.nvrt,nNodes) )
    if kbp2 is None:
        kbp2 = np.zeros( (nNodes,),dtype=int )
    ztmp[:,:] = np.nan

    iwet = eta+dp >= self.h0

    nz = max( self.nz, 1 ) # workaround for cases with no z levels
    hc = self.h_c # surf layer depth
    hs = self.h_s # s-z transition depth
    hmod = np.minimum( dp, hs )
    ztmp[nz-1,iwet] = -hmod[iwet]
    ztmp[self.nvrt-1,iwet] = eta[iwet]

    # check validity of v.grid (positive surf. layer depth)
    nSigma = self.nvrt-self.nz
    #etaMin = -hc-(hmod-hc)*self.cs[nSigma-1]/self.sigma[nSigma-1]
    #idx = np.logical_and( eta <= etaMin, iwet )
    #if idx.any() :
      ## raise error
      #errStr = ' '.join( ['Error in v.coords: choose a larger h_c :',
                          #str(hc), str(eta[idx][0]),str(etaMin[idx][0])] )
      #raise Exception( errStr )

    # case 1: dp <= hc, shallower than fine surface layer
    idx1 = np.logical_and( hmod <= hc, iwet )
    H1 = (hmod[idx1]+eta[idx1])
    eta1 = eta[idx1]
    # case 2: deeper than surf layer (normal case)
    idx2 = np.logical_and( hmod > hc, iwet )
    eta2 = eta[idx2]
    hmod2 = hmod[idx2]-hc
    for k in xrange(nz, self.nvrt-1) :
      kin = k - nz +1
      sigma = self.sigma[kin]
      ztmp[k,idx1] = sigma * H1 + eta1 # NOTE slow
      ztmp[k,idx2] = eta2*(1+sigma) + hc*sigma + hmod2*self.cs[kin] # NOTE slow
    # find bottom indices
    # pure sigma part
    idx = dp <= hs
    kbp2[idx] = nz-1
    # find z-coordinates
    idx_zpart = (dp > hs)
    for k in xrange(nz-1) :
      idx = idx_zpart & ( -dp >= self.ztot[k] ) & ( -dp < self.ztot[k+1] )
      kbp2[idx] = k
      ztmp[k,idx] = -dp[idx]
      idx = idx_zpart & ( -dp < self.ztot[k] )
      ztmp[k,idx] = self.ztot[k]
    ## extend z coordinates for shaved cells
    #hGapShaved = 0.5
    #for k in range(nz,-1,-1) :
      #idx = kbp2 == k
      #if idx.any() :
        #for kk in range(k-1,-1,-1) :
          #ztmp[kk,idx] = ztmp[kk+1,idx] - hGapShaved
    ztmp = np.ma.masked_invalid(ztmp)
    #print '  done.'
    return ztmp, kbp2, iwet


"""
lopezj - 08/06/2012
"""
class Object(object):
    pass

def readHGrid(path, gridonly=False):
    """ 
        Reads in an hgrid.gr3 file and places info into an object 
        
        grid - Object with grid information
    """
    try:
        f = open(path,"r")

        # Read the header information   
        header = f.readline()       
        [nElems, nNodes] = f.readline().rstrip("\n").split()
        nElems = int(nElems); nNodes = int(nNodes)

        # Get all the node information
        nodes = np.zeros((nNodes,4))
        for node in range(nNodes):
            [n, x, y, z] = f.readline().rstrip("\n").split()
            nodes[node,0] = int(n)+1
            nodes[node,1] = float(x)
            nodes[node,2] = float(y)
            nodes[node,3] = float(z)

        # Get the connectivity table
        elems = np.zeros((nElems,4))
        for elem in range(nElems):
            [n, x, a, b, c] = f.readline().rstrip("\n").split()
            elems[elem,0] = int(n)+1
            elems[elem,1] = int(a)
            elems[elem,2] = int(b)
            elems[elem,3] = int(c)

        # Get boundary node information
        if gridonly == False:
            obNodes = _readBndryNodes(f) # Open boundaries
            lbNodes = _readBndryNodes(f) # Land boundaries
            ibNodes = _readBndryNodes(f) # Island boundaries

        # Stuff everything into a generic object 
        grid = Object()
        grid.nodes = nodes
        grid.elems = elems
        if gridonly == False:
            grid.obNodes = obNodes
            grid.lbNodes = lbNodes
            grid.ibNodes = ibNodes
        return grid

    except IOError:
        raise
        
# get the grid  and startdate
argparser = argparse.ArgumentParser()
argparser.add_argument('hgrid', help='Name of hgrid.gr3 file')
argparser.add_argument('vgrid', help='Name of vgrid.in file')
argparser.add_argument('depthfile', help='Name of ocean_depths.nc file')
argparser.add_argument('basedir', help='Base directory for LiveOcean files')
argparser.add_argument('outdir', help='Output directory for CMOP nudging files')
argparser.add_argument('startdate', help='YYYY-MM-DD start date of nudging file')
argparser.add_argument('--ndays', type=int, help='number of days of forecast, defaults to 3',default=3)
args = argparser.parse_args()

hgrid = args.hgrid
vgrid = args.vgrid
depthfile = args.depthfile
startdate = args.startdate
basedir = args.basedir
outdir = args.outdir
ndays = args.ndays

try:
    grid = readHGrid(hgrid, True)
    x = grid.nodes[:,1]
    y = grid.nodes[:,2]
    d = grid.nodes[:,3]
    nnodes = len(x)
except:
    print('''Unable to open file %s''' % (hgrid,))
    sys.exit(1)
    
try:
    vCoords = vc.fromVGridFile(vgrid)
except Exception as e:
    print(e)
    
try:
    # construct vertical grid
    print('Vert grid')
    Z, kbp2, iwet = vCoords.computeVerticalCoordinates(d * 0.0,d)
except Exception as e:
    print(e)
    print('''Unable to open vertical grid file %s''' % (vgrid,))
    sys.exit(1)

nvrt = vCoords.nvrt

print(hgrid, startdate, nnodes)

print('ocean depths')
pncid = nc.Dataset(depthfile, 'r');
xrho = pncid.variables['lon_rho'][:]
yrho = pncid.variables['lat_rho'][:]

# PM Edit
#depths = pncid.variables['depths'][:]
depths = pncid.variables['h'][:]

xr=xrho[:]
yr=yrho[:]

depn = depths.shape[0]

# ROMS depths interpolated to SELFE grid
dep=np.zeros((depn,nnodes))
for i in range(0,depn):
    fz=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],depths[i,:,:].T,kx=1,ky=1)
    dep[i,:]=fz.ev(x,y)

(year, month, day) = startdate.split('-')
year = int(year)
month = int(month)
day = int(day)
t1 = parser.parse(startdate+'T00:00:00-00')
t2 = parser.parse(startdate+'T00:00:00-08')
dt= t2 -t1
startf = dt.total_seconds() / 3600.0
nfiles = ndays * 24 + 1
basedir += '''f%4d.%02d.%02d''' % (year,month,day)
print(basedir)
outname = '%s/cmop_%%s_nu.%4d-%02d-%02d.in' % (outdir,year,month,day)
sf = FortranFile(outname % 'salt', 'w')
tf = FortranFile(outname % 'temp', 'w')
temp_prev=np.zeros((depn,nnodes))
salt_prev=np.zeros((depn,nnodes))
print('processing')
for i in range(nfiles):
#for i in range(3):
# there is no file #1 use #2 twice
    tfile = startf + i + 1 #add 1 because there is an offset between file # and time step
    ncfile = '''%s/ocean_his_%04d.nc''' % (basedir, tfile,)
    print(ncfile)
    otime = i * 3600
    #for i in range(len(x)):
    temp_new=np.zeros((depn,nnodes))
    salt_new=np.zeros((depn,nnodes))
    try:
        psncid = nc.Dataset(ncfile, 'r');
        ocean_time = psncid.variables['ocean_time'][0]
        salt = psncid.variables['salt'][0]
        temp = psncid.variables['temp'][0]
        psncid.close()
        for j in range(0,depn):
            ft=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],temp[j,:,:].T,kx=1,ky=1)
            temp_new[j,:]=ft.ev(x,y)
            fs=interpolate.RectBivariateSpline(xrho[0,:],yrho[:,0],salt[j,:,:].T,kx=1,ky=1)
            salt_new[j,:]=fs.ev(x,y)
        temp_new[np.where(salt_new>35)]=float('nan')
        salt_new[np.where(salt_new>35)]=float('nan')
        good=~np.isnan(salt_new[0,:])
        bad=np.isnan(salt_new[0,:])
        xygood=np.array((x[good],y[good])).T
        xybad=np.array((x[bad],y[bad])).T    
        ind=np.arange(nnodes)
        map_tr=ind[good]
        kdt=cKDTree(zip(x[good],y[good]))
        dists,neighs=kdt.query(xybad)
        min_index=map_tr[neighs]
        ix_min=np.unravel_index(min_index,x.shape)
        salt_new[:,bad]=salt_new[:,ix_min][:,0,:]
        temp_new[:,bad]=temp_new[:,ix_min][:,0,:]
        salt_prev = salt_new
        temp_prev = temp_new
    except:
        print('Exception in loading new, using previous time step')
        salt_new = salt_prev
        temp_new = temp_prev
    print(i+1, ' of ', nfiles, tfile, otime, ocean_time, ncfile)
    sf.write_record(np.array([otime], dtype=np.float32))
    tf.write_record(np.array([otime], dtype=np.float32))
    for j in range(len(x)):
        s1=np.zeros(nvrt)
        t1=np.zeros(nvrt)
        ss=interp(Z[:,j],dep[:,j],salt_new[:,j])
        tt=interp(Z[:,j],dep[:,j],temp_new[:,j])
        s1[Z[:,j].mask==False]=ss[Z[:,j].mask==False]
        t1[Z[:,j].mask==False]=tt[Z[:,j].mask==False]
        sf.write_record(s1.astype(np.float32))
        tf.write_record(t1.astype(np.float32))

import os
import sys

# add .. to PYTHONPATH
path = os.path.realpath(__file__)
#print('My path', path)
basedir = os.path.dirname(os.path.dirname(path))

#print('PYTHONPATH is', os.environ['PYTHONPATH'])

#print('Adding basedir', basedir, 'to PYTHONPATH')
sys.path.append(basedir)

# add ../blind and ../util to PATH
os.environ['PATH'] += ':' + os.path.join(basedir, 'blind')
os.environ['PATH'] += ':' + os.path.join(basedir, 'util')

from astropy.io import fits as pyfits

import tempfile
import traceback
import multiprocessing
import time
import re
import gzip
import zipfile
import math
import time

import logging

#from astrometry.util import image2pnm
from astrometry.util.filetype import filetype_short
from astrometry.util.run_command import run_command

from astrometry.util.util import Tan
from astrometry.util import util as anutil

import argparse

SCALEUNITS_CHOICES = (
    ('arcsecperpix', 'arcseconds per pixel'),
    ('arcminwidth' , 'width of the field (in arcminutes)'), 
    ('degwidth' , 'width of the field (in degrees)'),
    ('focalmm'     , 'focal length of the lens (for 35mm film equivalent sensor)'),
)
scaleunits_default = 'degwidth'

SCALETYPE_CHOICES = (
    ('bounds', 'bounds'),
    ('estimate', 'estimate +/- error'),
)

wcsfile = 'wcs.fits'
corrfile = 'corr.fits'
axyflags = []

def desc_from_keyvalues(kvlist):
    return ', '.join([': '.join(x) for x in kvlist])

parser = argparse.ArgumentParser(description='Simple Integrated pipline')

parser.add_argument('--scaleunits', choices=SCALEUNITS_CHOICES, default=scaleunits_default, help=desc_from_keyvalues(SCALEUNITS_CHOICES))
parser.add_argument('--scaletype', choices=[x[0] for x in SCALETYPE_CHOICES], default=SCALETYPE_CHOICES[1][0], help=desc_from_keyvalues(SCALETYPE_CHOICES))
parser.add_argument('--parity', choices = ['pos', 'neg', 'both'], default='both', help='consider one parity case or both')
parser.add_argument('--scale-lower', type=float, default=0.1, help='lower bound of scale, for bounds scaletype model')
parser.add_argument('--scale-upper', type=float, default=180, help='upper bound of scale, for bounds scaletype model')
parser.add_argument('--scale-estimate', type=float, default=None, help='estimated scale, for estimated scaletype model')
parser.add_argument('--scale-error', type=float, default=None, help='estimated scale error, for estimated scaletype model')
parser.add_argument('--pixel-error', type=float, default=None, help='pixel error')
parser.add_argument('--ra', type=float, default=None, help='estimated RA')
parser.add_argument('--dec', type=float, default=None, help='estimated DEC')
parser.add_argument('--search-radius', type=float, default=5, help='search radius, in degrees')
parser.add_argument('--tweak-order', type=int, default=2, help='tweak order...? help needs improvement')
parser.add_argument('--downsample', type=int, default=2, help='allow downsampling before solve by this factor')
# NOTE, these are ONLY to hold user-set (via API) image width/height;
# they OVERRIDE the actual size of the image.  ONLY valid for xylists.
parser.add_argument('--image-width', type=int, default=0, help='image width for override, in xylist input (without explicit-image)')
parser.add_argument('--image-height', type=int, default=0, help='image height for override, in xylist input (without explicit-image)')
parser.add_argument('--use-sextractor', action='store_true')
parser.add_argument('--crpix-center', action='store_true')
parser.add_argument('--invert', action='store_true')
parser.add_argument('-v', '--verbosity', choices=['debug', 'info', 'warning', 'error', 'critical'], default='info', help='log verbosity level')
parser.add_argument('inputfile', type=str, help='input file (fits or xylist)')

args = parser.parse_args()

logging.basicConfig(filename='output.log',level=getattr(logging, args.verbosity.upper()))
logger = logging.getLogger()
logger.setLevel(getattr(logging, args.verbosity.upper()))

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(getattr(logging, args.verbosity.upper()))
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

fh = logging.FileHandler('{}.log'.format(os.path.splitext(os.path.basename(args.inputfile))[0]))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

if args.scaletype == 'estimated' and (args.scale_estimate is None or args.scale_error is None):
    parser.error('scaletype=ev requires --scale-estimate=estimated_scale --scale-error=estimated_error')

scale_lower = args.scale_lower
scale_upper = args.scale_upper
scale_estimiate = args.scale_estimate
scale_error = args.scale_error
image_width = args.image_width
image_height = args.image_height

slo, shi = scale_lower, scale_upper
if args.scaletype == 'ev':
    slo,shi = (scale_estimate * (1.0 - scale_error / 100.0), scale_estimate * (1.0 + scale_error / 100.0))

jobdir = '.'
axyfn = 'result.axy'
axypath = os.path.join(jobdir, axyfn)

axyargs = {
    '--out': axypath,
    '--scale-low': slo,
    '--scale-high': shi,
    '--scale-units': args.scaleunits,
    '--wcs': wcsfile,
    '--corr': corrfile,
    '--rdls': 'rdls.fits',
    '--pixel-error': args.pixel_error,
    '--downsample': args.downsample,
    # tuning-up maybe fixed; if not, turn it off with:
    #'--odds-to-tune': 1e9,

    # Other things we might want include...
    # --invert
    # -g / --guess-scale: try to guess the image scale from the FITS headers
    # --crpix-x <pix>: set the WCS reference point to the given position
    # --crpix-y <pix>: set the WCS reference point to the given position
    # -w / --width <pixels>: specify the field width
    # -e / --height <pixels>: specify the field height
    # -X / --x-column <column-name>: the FITS column name
    # -Y / --y-column <column-name>
}

if args.ra and args.dec and args.search_radius:
    axyargs['--ra'] = args.ra
    axyargs['--dec'] = args.dec
    axyargs['--radius'] = args.radius

from pipes import quote

image_height = args.image_height
image_width = args.image_width
if False: #hasattr(img,'sourcelist'):
    # image is a source list; use --xylist
    axyargs['--xylist'] = quote( args.inputfile )
    w,h = img.width, img.height
    if args.image_width:
        w = args.image_width
    if args.image_height:
        h = args.image_height
    axyargs['--width' ] = w
    axyargs['--height'] = h
else:
    axyargs['--image'] = quote( args.inputfile )
    with pyfits.open(args.inputfile) as img:
        image_height, image_width = img[0].shape

if args.parity != 'both':
    axyargs['--parity'] = args.parity

axyflags = []
if args.tweak_order == 0:
    axyflags.append('--no-tweak')
else:
    axyargs['--tweak-order'] = '%i' % args.tweak_order

if args.use_sextractor:
    axyflags.append('--use-sextractor')

if args.crpix_center:
    axyflags.append('--crpix-center')

if args.invert:
    axyflags.append('--invert')
    
cmd = 'augment-xylist '

def join(terms):
    return ' '.join(terms)

cmd = 'augment-xylist {}'.format(join([join([' '.join([str(y) for y in x]) for x in axyargs.iteritems()]), join(axyflags)]))
logger.info('running: ' + cmd)
xylist_start = time.time()
(rtn, out, err) = run_command(cmd, tee=True)
xylist_end = time.time()
if rtn:
    logger.info('out: ' + out)
    logger.info('err: ' + err)
    logger.critical('augment-xylist failed: rtn val: {} err: {}'.format(rtn, err))
    sys.exit(-1)
logger.info('created axy file: {} in {:0.2f} seconds'.format(axypath, xylist_end-xylist_start))

cmd = 'astrometry-engine {}'.format(axypath)
astrometry_engine_start = time.time()
(rtn, out, err) = run_command( cmd, tee=True )
astrometry_engine_end = time.time()
if rtn:
    logger.info('astrometry-engine output: ' + out)
    logger.info('astrometry-engine err: ' + err)
    logger.critical('astrometry-engine failed: rtn val {} err: {}'.format(rtn, err))
    sys.exit(-1)

wcsfn = os.path.join(jobdir, wcsfile)
logger.debug('Checking for WCS file {}:'.format(wcsfn))
if os.path.exists(wcsfn):
    logger.debug('WCS file exists, astrometry-engine solved in {:0.2f} seconds'.format(astrometry_engine_end - astrometry_engine_start))
logger.info('astrometry-engine solved in {:0.2f} seconds'.format(astrometry_engine_end - astrometry_engine_start))

# Parse the wcs.fits file
wcs = Tan(wcsfn, 0)
model_solve_text = ('<TanWCS: CRVAL (%f, %f)' % (wcs.crval[0], wcs.crval[1]) +
        ' CRPIX (%f, %f)' % (wcs.crpix[0], wcs.crpix[1]) +
        ' CD (%f, %f; %f %f)' % tuple(wcs.cd) +
        ' Image size (%d, %d)>' % (image_width, image_height)
        )
logger.info('Created Tan-WCS solution!: {}'.format( model_solve_text ))

# Find field's healpix nside and index

ra,dec = wcs.radec_center()
radius = wcs.pixel_scale() * math.hypot(wcs.imagew, wcs.imageh)/2.0 / 3600.0

nside = anutil.healpix_nside_for_side_length_arcmin(radius*60)
nside = int(2**round(math.log(nside, 2)))
nside = max(1, nside)
healpix = anutil.radecdegtohealpix(ra, dec, nside)

# Find bounds for the Calibration object.
r0,r1,d0,d1 = wcs.radec_bounds()
# Find cartesian coordinates
ra *= math.pi/180
dec *= math.pi/180
tempr = math.cos(dec)
x = tempr*math.cos(ra)
y = tempr*math.sin(ra)
z = math.sin(dec)
r = radius/180*math.pi

output="""
Solution:
    ra: {ra}
    dec: {dec}
    x: {x} y: {y} z: {z}
    r: {r}
    RA bounds: [{ralo}, {rahi}]
    DEC bounds: [{declo}, {dechi}]
""".format(ra=ra, dec=dec, x=x, y=y, z=z, r=r, ralo=r0, rahi=r1, declo=d0, dechi=d1)

print(output)
import os
import sys

##
## TLSMD application settings
##

## program info
VERSION       = "0.9.0 Personal Edition"
RELEASE_DATE  = "1 Jul 2007"
AUTHOR        = "Jay Painter"
EMAIL         = "jay.painter@gmail.com"


## Root directory of TLSMD installation
TLSMD_ROOT = os.getcwd()

## Path to GNU Plot executeable
GNUPLOT_PATH = '/usr/bin/gnuplot'


## BEGIN: CONFIGURATION PATHS AND URLS
TLSMD_BASE_URL         = "/tlsmd"
WEBTLSMDD              = "http://localhost:10100"
WEBTLSMDD_DATABASE     = os.path.join(TLSMD_ROOT, "run", "webtlsmd.db")
ADMIN_PASSWORD_FILE    = os.path.join(TLSMD_ROOT, "run", "admin-password")
TRACEBACK_EMAIL        = "tlsmdtraceback"
LOG_PATH               = os.path.join(TLSMD_ROOT, "run", "tlsmd_runlog.txt")
## END: CONFIGURATION PATHS AND URLS

## derived paths
TLSMD_PROGRAM_PATH     = os.path.join(TLSMD_ROOT, "bin", "tlsmd.py")
GNUPLOT_FONT           = os.path.join(TLSMD_ROOT, "fonts/LucidaSansOblique.ttf")
REFINEPREP_URL         = "/tlsmd/refineprep.html"
TLSMD_WORK_DIR         = os.path.join(TLSMD_ROOT, "static", "jobs")
TLSMD_WORK_URL         = "/tlsmd/static/jobs"
JMOL_DIR               = "/static/java"

## maximum number of residues in a chain
MAX_CHAIN_LENGTH = 1700

## the isoprobability contour level for all visualizations
ADP_PROB = 50

## number of TLS partitons for each chain
NPARTS = 20

## the pixel width of the TLS visualization rendered ray traces
VIS_WIDTH = 640

## the JMol viewer is a square window, generated with this pixel size
JMOL_SIZE = 600


##
## DJANGO Settings Below
##

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    ('Jay Painter', 'jay.painter@gmail.com'),
)


EMAIL_SUBJECT_PREFIX = 'TLSMD Server'
SERVER_EMAIL = 'jay.painter@gmail.com'


MANAGERS = ADMINS

## DATABASE CONFIGURATION FOR SQLITE3
## DATABASE_ENGINE = 'sqlite3' # 'postgresql', 'mysql', 'sqlite3' or 'ado_mssql'.
## DATABASE_NAME = os.path.join(TLSMD_ROOT, 'run', 'django-webtlsmd.db') # Or path to database file if using sqlite3.
## DATABASE_USER = ''             # Not used with sqlite3.
## DATABASE_PASSWORD = ''         # Not used with sqlite3.
## DATABASE_HOST = ''             # Set to empty string for localhost. Not used with sqlite3.
## DATABASE_PORT = ''             # Set to empty string for default. Not used with sqlite3.

## DATABASE CONFIG FOR MYSQL
DATABASE_ENGINE = 'mysql'      # 'postgresql', 'mysql', 'sqlite3' or 'ado_mssql'.
DATABASE_NAME = 'tlsmd'        # Or path to database file if using sqlite3.
DATABASE_USER = 'root'         # Not used with sqlite3.
DATABASE_PASSWORD = ''         # Not used with sqlite3.
DATABASE_HOST = ''             # Set to empty string for localhost.
DATABASE_PORT = ''             # Set to empty string for default


# Local time zone for this installation. All choices can be found here:
# http://www.postgresql.org/docs/current/static/datetime-keywords.html#DATETIME-TIMEZONE-SET-TABLE
TIME_ZONE = 'America/Chicago'

# Language code for this installation. All choices can be found here:
# http://www.w3.org/TR/REC-html40/struct/dirlang.html#langcodes
# http://blogs.law.harvard.edu/tech/stories/storyReader$15
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# Absolute path to the directory that holds media.
# Example: "/home/media/media.lawrence.com/"
MEDIA_ROOT = ''

# URL that handles the media served from MEDIA_ROOT.
# Example: "http://media.lawrence.com"
MEDIA_URL = ''

# URL prefix for admin media -- CSS, JavaScript and images. Make sure to use a
# trailing slash.
# Examples: "http://foo.com/media/", "/media/".
ADMIN_MEDIA_PREFIX = '/media/'

# Make this unique, and don't share it with anybody.
SECRET_KEY = '5_*ndl@fd1j&z=plw^7hlo-fdkx#44@(^(wqvo*m4@r#ofya9l'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.load_template_source',
    'django.template.loaders.app_directories.load_template_source',
#     'django.template.loaders.eggs.load_template_source',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.middleware.doc.XViewMiddleware',
)

ROOT_URLCONF = 'tlsmdlib.urls'

TEMPLATE_DIRS = (
    os.path.join(TLSMD_ROOT, 'templates'),
)

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'tlsmdlib',
)

"""
pyXMM - Python interface to XMM SAS for data reduction and extraction of
        XMM-Newton products

v2.0 - 15/03/2017 - D.R. Wilkins

Data reduction/extraction for a single OBSID is handled by EPICExtractor class
"""
import os
import glob
import re
import sys
import subprocess
import math
import pyfits

from .epic_extractor import *
from .gti import *
from .util import *
from .spec_util import *
from .spec_combine import *

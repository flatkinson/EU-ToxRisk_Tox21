import warnings
import logging
import os
import sys
import re
import json
from math import ceil
import numpy as np
import pandas as pd
import requests
logging.getLogger('requests').setLevel('WARNING')
from ipywidgets import HTML

from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = (450, 200)
PandasTools.RenderImagesInAllDataFrames()

import make_logger

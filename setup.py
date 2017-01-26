import warnings
import logging
import sys
import re
import json
import csv
from math import ceil
import pandas as pd
import numpy as np
import requests
logging.getLogger('requests').setLevel('WARNING')
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = (450, 200)
PandasTools.RenderImagesInAllDataFrames()
from ipywidgets import HTML

import make_logger

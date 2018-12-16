#!/usr/bin/env python3

import os
import numpy as np

Vgrid = np.linspace(0,1,11)

for Vg in Vgrid:
    os.system('python3 Step1_code.py %1.1f' %Vg)


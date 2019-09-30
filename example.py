import os
import sys
# cwd = os.getcwd()
# if(sys.version_info[0]==3 and sys.version_info[1]==5):
#     sys.path.insert(0, os.path.normpath(os.path.join(cwd,'lib','windows','cp35')))
# elif(sys.version_info[0]==3 and sys.version_info[1]==6):
#     sys.path.insert(0, os.path.normpath(os.path.join(cwd,'lib','windows','cp36')))
# elif(sys.version_info[0]==3 and sys.version_info[1]==7):
#     sys.path.insert(0, os.path.normpath(os.path.join(cwd,'lib','windows','cp37')))

from lsreader import D3plotReader
from lsreader import DataType as dt

# father_path = os.path.abspath(os.path.dirname(cwd)+os.path.sep+".")
# data_path = os.path.join(father_path, 'example_data', 'd3plot')
# data_path = os.path.abspath(data_path)
# dr = D3plotReader(data_path)

# shell_stress = dr.get_data(dt.D3P_SHELL_STRESS, ist = 0, ipt = 1)
# print(shell_stress[0].x())

# shell_eps = dr.get_data(dt.D3P_SHELL_EFFECTIVE_PLASTIC_STRAIN, ist = 0, ipt = 1)
# print(shell_eps[0])

# thickness = dr.get_data(dt.D3P_SHELL_THICKNESS, ist = 0, ipt = 1)
# print(thickness[0])

# solid_element = dr.get_data(dt.D3P_NUM_SOLID, ist = 0, ipt = 1)
# print(solid_element)


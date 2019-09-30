import struct
import array
import numpy as np


def read_stl(fname):
    f = open(fname, 'rb')
    rez = {'title': struct.unpack('80c', f.read(80))[0]}
    rez['polygons_count'] = struct.unpack('I', f.read(4))[0]
    rez['polygons'] = []
    for i in range(rez['polygons_count']):
        polygon = struct.unpack('12f', f.read(48))
        pol = {'normal': polygon[:3],
               'points': [polygon[3+i*3: 3+(i+1)*3] for i in range(3)],
               'flag': struct.unpack('H', f.read(2))[0]
               }
        rez['polygons'].append(pol)
    f.close()
    return rez


model = read_stl('1.stl')
print(model['polygons_count'])

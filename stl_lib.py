import struct
import array
import numpy as np
from math import sqrt
import json
from os.path import exists


class Point(object):
    def __init__(self, x, y, z):
        """
        Point(x,y,z) - create point(vector) in 3D
        """
        self.x = x
        self.y = y
        self.z = z

    def norm(self):
        """
        Length of the vector
        """
        return sqrt(self.x*self.x+self.y*self.y+self.z*self.z)

    def normalized(self):
        """
        Return vector of unity length
        """
        l = self.norm()
        if l == 0:
            return self
        p = Point(self.x/l, self.y/l, self.z/l)
        return p

    def __repr__(self):
        return '(%f,%f,%f)' % (self.x, self.y, self.z)

    def dot(self, p2):
        rez = self.x*p2.x+self.y*p2.y+self.z*p2.z
        return rez

    def cross(self, p2):
        x = self.y*p2.z-self.z*p2.y
        y = self.z*p2.x-self.x*p2.z
        z = self.x*p2.y-self.y*p2.x
        return Point(x, y, z)

    def scaled(self, a):
        return Point(self.x*a, self.y*a, self.z*a)

    def __add__(self, p2):
        return Point(self.x+p2.x, self.y+p2.y, self.z+p2.z)

    def __sub__(self, p2):
        return Point(self.x-p2.x, self.y-p2.y, self.z-p2.z)

    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)

    def __rmul__(self, r):
        return self.scaled(r)

    def __mul__(self, r):
        return self.scaled(r)

    def shuff(self, p1, p2):
        return (self.cross(p1)).dot(p2)

    def _getData(self):
        return np.array([self.x, self.y, self.z])

    def _setData(self, coords):
        coords = list(coords) + [0]*(3-len(coords))
        self.x, self.y, self.z = coords

    data = property(_getData, _setData)

    def angle(self, p2, direction_vector=None):
        p1_u = self.normalized()
        p2_u = p2.normalized()
        rez = np.arccos(
            np.clip(np.dot(p1_u.data, p2_u.data), -1.0, 1.0))*180./np.pi
        if direction_vector:
            if self.cross(p2).dot(direction_vector) < 0:
                rez = 360-rez
        return rez


def line_intersect_plane(line_p0, line_p1, plane_point, plane_normal, full_line=False):
    """
    function find the intersection between line, defined by points line_p0 and line_p1
    and plane, defined by point plane_point and normal plane_normal.
    if line is parallel to plane and there is no intersection function returns 0.
    if full_line = True, line is extended to intersection,
    else point returns only if segment p0-p1 intersects plane.
    """
    v1 = (line_p1-line_p0).dot(plane_normal)
    if v1 == 0:
        return None
    v2 = (plane_point-line_p0).dot(plane_normal)
    rI = v2/v1
    lp = line_p0+rI*(line_p1-line_p0)
    if full_line:
        return lp
    else:
        if 0 <= rI <= 1:
            return lp
        else:
            return None


class Polygon(object):
    """triangular poligon"""

    def __init__(self, points, normal, flag):
        super(Polygon, self).__init__()
        self.points = points
        self.normal = normal
        self.flag = flag

    def __repr__(self):
        rez = 'Normal:' + repr(self.normal)+'\n'
        rez += 'Points:\n'
        for i in range(3):
            rez += repr(self.points[i])+'\n'
        rez += 'flag:' + str(self.flag)+'\n'
        return rez

    def bbox(self):
        xs = [self.points[i].x for i in range(3)]
        ys = [self.points[i].y for i in range(3)]
        zs = [self.points[i].z for i in range(3)]
        return [(min(xs), max(xs)),
                (min(ys), max(ys)),
                (min(zs), max(zs))
                ]

    def center(self):
        return np.sum(self.points)*0.3333333333333

    def plane_intersection(self, plane_point, plane_normal):
        def plane(p): return plane_normal.dot(p-plane_point)
        rez = np.sign([plane(pp) for pp in self.points])
        if np.all(rez == 0):
            return lambda r: self.center()
        if rez[0] == rez[1] == rez[2]:
            return None
        _, idx, counts = np.unique(rez, return_index=True, return_counts=True)
        min_idx = idx[counts.argmin()]
        p0 = self.points[min_idx]
        if rez[min_idx] == 0:
            return p0
        idxs = [0, 1, 2]
        idxs.pop(min_idx)
        p1 = self.points[idxs[0]]
        p2 = self.points[idxs[1]]
        pstart = line_intersect_plane(p0, p1, plane_point, plane_normal)
        pend = line_intersect_plane(p0, p2, plane_point, plane_normal)
        return (pstart+pend)*0.5


class StlModel(object):
    def __init__(self, fname):
        print('reading stl model from '+fname)
        self.fname = fname
        self.polygons = read_stl(fname)['polygons']
        print('done reading slt model')
        self.boundbox = None

    def devide_on_blocks(self, dir=0, minc=0, maxc=1, ndivs=2):
        z = np.linspace(minc, maxc, ndivs+1)
        self.blocks = {}
        for i in range(1, len(z)):
            block_name = '{} {}'.format(z[i-1], z[i])
            self.blocks[block_name] = []
            for j, p in enumerate(self.polygons):
                if z[i-1] <= p.center().data[dir] < z[i]:
                    self.blocks[block_name].append(j)

    def bbox(self):
        if self.boundbox:
            return self.boundbox
        minc = [1e10, 1e10, 1e10]
        maxc = [-1e10, -1e10, -1e10]
        for pol in self.polygons:
            for p in pol.points:
                for i in range(3):
                    minc[i] = min(minc[i], p.data[i])
                    maxc[i] = max(maxc[i], p.data[i])
        self.boundbox = list(zip(minc, maxc))
        return self.boundbox

    def save_blocks(self):
        if self.blocks:
            json.dump(self.blocks, open(self.fname+'.blks', 'w'))

    def read_blocks(self):
        if exists(self.fname+'.blks'):
            self.blocks = json.load(open(self.fname+'.blks', 'r'))


def read_stl(fname):
    f = open(fname, 'rb')
    rez = {'title': struct.unpack('80c', f.read(80))[0]}
    rez['polygons_count'] = struct.unpack('I', f.read(4))[0]
    rez['polygons'] = []
    for i in range(rez['polygons_count']):
        polygon = struct.unpack('12f', f.read(48))
        pol = Polygon(normal=Point(*polygon[:3]),
                      points=[Point(*polygon[3+i*3: 3+(i+1)*3])
                              for i in range(3)],
                      flag=struct.unpack('H', f.read(2))[0]
                      )
        rez['polygons'].append(pol)
    f.close()
    return rez


if __name__ == '__main__':
    pass

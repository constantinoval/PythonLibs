import struct
import array
import numpy as np


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

    def __lmul__(self, r):
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


def line_intersect_plane(line_p0, line_p1, plane_point, plane_normal):
    v1 = (line_p1-line_p0).dot(plane_normal)
    if v1 == 0:
        return 0
    v2 = (plane_point-line_p0).dot(plane_normal)
    rI = v2/v1
    return line_p0+rI*(line_p1-line_p0)


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
        return rez


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
    # model = read_stl('1.stl')
    # print(model['polygons'])
    print(line_intersect_plane(Point(0, 1, 0),
                               Point(1, 1, 0),
                               Point(2, 0, 0),
                               Point(1, 0, 0)))

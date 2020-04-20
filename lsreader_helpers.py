import lsreader as ls
import numpy as np


class d3plot_data(object):
    def __init__(self, path_to_d3plot):
        self.reader = ls.D3plotReader(path_to_d3plot)
        self.parts = None
        self.solids = None
        self.nodes = None
        self.shells = None

    def _get_parts(self):
        dtype = {'names': ['index', 'ID', 'name'],
                 'formats': [np.int32, np.int32, 'O']}
        nparts = self.reader.get_data(ls.DataType.D3P_NUM_PARTS)
        IDs = self.reader.get_data(ls.DataType.D3P_PART_IDS)
        names = [str(self.reader.get_data(ls.DataType.D3P_PART_NAME, ipart=i)).strip()
                 for i in range(nparts)]
        return np.array(
            list(zip(list(np.arange(nparts)), list(IDs), names)), dtype=dtype)

    def read_parts(self):
        self.parts = self._get_parts()

    def get_number_of_steps(self):
        return self.reader.get_data(ls.DataType.D3P_NUM_STATES)

    def get_times(self):
        return list(self.reader.get_data(ls.DataType.D3P_TIMES))

    def get_part_velocity(self, partID):
        if not self.parts:
            self.read_parts()
        if not partID in self.parts['ID']:
            print(f'Part {partID} not found')
            return None
        pi = self.parts[self.parts['ID'] == partID]['index'][0]
        v = {'x': [], 'y': [], 'z': []}
        for s in range(self.get_number_of_steps()):
            _ = self.reader.get_data(
                ls.DataType.D3P_PART_VELOCITY, ist=s, ipart=int(pi))
            v['x'].append(_.x())
            v['y'].append(_.y())
            v['z'].append(_.z())
        return v

    def mesh_info(self):
        self.nSolids = self.reader.get_data(ls.DataType.D3P_NUM_SOLID)
        self.nSPH = self.reader.get_data(ls.DataType.D3P_NUM_SPH)
        self.nShells = self.reader.get_data(ls.DataType.D3P_NUM_SHELL)
        self.nTShells = self.reader.get_data(ls.DataType.D3P_NUM_TSHELL)
        self.nNodes = self.reader.get_data(ls.DataType.D3P_NUM_NODES)
        self.nParts = self.reader.get_data(ls.DataType.D3P_NUM_PARTS)

    def read_solids(self):
        dtype = {'names': ['index', 'ID', 'mat', 'n1', 'n2', 'n3', 'n4', 'n5',
                           'n6', 'n7', 'n8', 'n9', 'n10'],
                 'formats': [np.int32]*13}
        nsolids = self.reader.get_data(ls.DataType.D3P_NUM_SOLID)
        IDs = self.reader.get_data(ls.DataType.D3P_SOLID_IDS)
        els = self.reader.get_data(ls.DataType.D3P_SOLID_CONNECTIVITY_MAT)
        mat = []
        n = [[], [], [], [], [], [], [], [], [], []]
        for e in els:
            mat.append(e.mat()-1)
            for i in range(10):
                n[i].append(e.node(i))
        self.solids = np.array(
            list(zip(list(np.arange(nsolids)), IDs, mat, *n)), dtype=dtype)

    def read_nodes(self):
        self.nNodes = self.reader.get_data(ls.DataType.D3P_NUM_NODES)
        dtype = {'names': ['index', 'ID'],
                 'formats': [np.int32, np.int32]}
        IDs = self.reader.get_data(ls.DataType.D3P_NODE_IDS)
        self.nodes = np.array(
            list(zip(list(np.arange(self.nNodes)), list(IDs))), dtype=dtype)

    def get_solids_by_partID(self, partID):
        if not self.parts:
            self.read_parts()
        if not self.solids:
            self.read_solids()
        return self.solids[self.solids['mat'] == self.parts[self.parts['ID'] == partID]['index'][0]]

    def get_nodes_by_patrID(self, partID):
        if not self.nodes:
            self.read_nodes()
        els = self.get_solids_by_partID(partID)
        ni = []
        for e in els:
            for n in ['n1', 'n2', 'n3', 'n4', 'n5',
                      'n6', 'n7', 'n8', 'n9', 'n10']:
                if e[n]:
                    ni.append(e[n])
        ni = list(set(ni))
        return self.nodes[ni]

    def get_node_coord(self, nID=-1, nIndex=-1, ist=0):
        if nID >= 0:
            if not self.nodes:
                self.read_nodes()
            nIndex = int(self.nodes[self.nodes['ID'] == nID]['index'][0])
        crds = self.reader.get_data(
            ls.DataType.D3P_NODE_COORDINATES, ist=ist)[nIndex]
        return [crds.x(), crds.y(), crds.z()]
        # return crds


if __name__ == '__main__':
    # import matplotlib.pylab as plt
    d = d3plot_data('d3plot')
    # v = d.get_part_velocity(20)
    # t = d.get_times()
    # plt.plot(t, v['x'])
    # plt.show()
    # d.read_nodes()
    #p20 = d.get_solids_by_partID(20)
    #n20 = d.get_nodes_by_patrID(20)
    print(d.get_node_coord(nIndex=116221))

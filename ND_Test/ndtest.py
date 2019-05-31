import numpy as np

class Simplex:
    def __init__(self):
        self.num_dimensions = 2
        self.points = []
        for d in range(self.num_dimensions + 1):
            p = []
            for dd in range(self.num_dimensions):
                p = p + [0.0]
            self.points = self.points + [p]
        self.points = np.array(self.points)

class Cell:
    def __init__(self):
        self.num_dimensions = 2
        self.points = []
        for d in range(2**self.num_dimensions):
            p = []
            for dd in range(self.num_dimensions):
                p = p + [0.0]
            self.points = self.points + [p]
        self.points = np.array(self.points)

class Grid:
    def __init__(self, _dimensions, _resolution):
        self.num_dimensions = len(_dimensions)
        self.dimensions = _dimensions
        self.resolution = _resolution
        self.points = self.generate_points([], self.dimensions, self.resolution)
        l = self.generate_cells([],self.resolution)
        self.cells = []

    def generate_points(self, total, dim, res):
        dh, *dt = dim
        rh, *rt = res

        if len(dt) == 0:
            ps = []
            for d in range(rh+1):
                ps = ps + [total + [d * (dh / rh)]]
            return ps


        ps = []
        for d in range(rh+1):
            ps = ps + [self.generate_points(total + [d*(dh/rh)], dt, rt)]
        return ps


    def int_to_bin_list(self, i):
        bstring = str(bin(i))
        return [int(a) for a in bstring[2:].zfill(self.num_dimensions)]

    def cell_index_list(self):
        return [self.int_to_bin_list(a) for a in range(2**self.num_dimensions)]

    def get_cell_indices(self, cell):
        return [(np.array(a) + np.array(cell)).tolist() for a in self.cell_index_list()]

    def generate_cells(self, cell, res):
        rh, *rt = res

        cells = []
        if len(rt) == 0:
            for r in range(rh):
                cells = cells + [self.get_cell_indices(cell + [r])]
            return cells

        for d in range(rh):
            cells = cells + [self.generate_cells(cell + [d], rt)]

        return cells


g = Grid([3.0,3.0,3.0], [3,3,3])

# ------------------------------------------------------------------------
# Copyright 2018,  Alberto Pietrini
# Statisticalhitfinder is distributed under the terms of the Simplified BSD License.
# -------------------------------------------------------------------------

import numpy as np


class Functions:
    def __init__(self, nr_workers, rank):
        self.size = nr_workers
        self.rank = rank

    def chunk_file(self, dimensions, flag):
        self.x, self.y, self.z = dimensions
        self.flag = flag

        # FRAMEWISE
        if self.flag == 'framewise':
            self.step = np.round(self.x/self.size)
            self.start = self.rank*self.step

            if self.rank != self.size - 1:
                self.end = self.step*(self.rank+1)
            else:
                self.end = self.x

            return self.start, self.end

        # PIXELWISE
        elif self.flag == 'pixelwise':
            self.step = np.round(self.y/self.size)
            self.start = self.rank*self.step

            if self.rank != self.size - 1:
                self.end = self.step*(self.rank+1)
                self.ext_dim = self.z*self.step
            else:
                self.end = self.y
                self.ext_dim = self.z*(self.end-self.start)

            return self.start, self.end, self.ext_dim

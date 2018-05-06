# ------------------------------------------------------------------------
# Copyright 2018,  Alberto Pietrini
# Statisticalhitfinder is distributed under the terms of the Simplified BSD License.
# -------------------------------------------------------------------------


# READ AND WRITE H5 FILES IN PARALLEL MODE

from mpi4py import MPI
import h5py

# mpi_def DEFINE ALL NEEDED TO CALL MPI

class mpi_h5rw:
    def __init__(self):
        self.comm = MPI.COMM_WORLD
        self.size = self.comm.size
        self.rank = self.comm.rank

    def createH5(self, filename, rw_flag):
        self.filename = filename
        self.rw_flag = rw_flag
        self.fileID = h5py.File(self.filename, self.rw_flag, driver='mpio', comm=self.comm)
        if rw_flag == 'w': self.fileID.atomic = True

        # WAIT FOR ALL WORKERS TO CREATE/READ ALL THE FILES
        self.comm.barrier()

        return self.fileID

    def create_ds(self, fileID, ds_name, ds_size, ds_type):
        # CREATE DATASET IN THE H5
        self.fileID = fileID
        self.ds_name = ds_name
        self.ds_size = ds_size
        self.ds_type = ds_type
        # ds_name HERE SHOULD BE GIVEN LIKE THAT: '/entry_1/data_1/data'
        self.ds = self.fileID.create_dataset(self.ds_name, self.ds_size, dtype=self.ds_type)

        # self.ds SHOULD BE USED TO PUT DATA IN, LIKE: self.ds[:] = data
        # OPTION self.ds.flush() SHOULD BE USED IN THE MAIN PROGRAM TO AVOID MEMORY ISSUES

        return self.ds

    def read_ds(self, fileID, ds_name):
        self.fileID = fileID
        self.ds_name = ds_name
        # ds_name HERE SHOULD BE GIVEN LIKE THAT: '/entry_1/data_1/data'
        self.ds = self.fileID[self.ds_name]

        return self.ds

    def ds_flush(self, fileID):
        self.fileID = fileID
        self.fileID.flush()
        #self.comm.barrier()

    def h5_close(self, fileID):
        self.fileID = fileID
        self.fileID.close()

    def mpi_params(self):
        return self.size, self.rank, self.comm

#!/usr/bin/env python
# ------------------------------------------------------------------------
# Copyright 2018,  Alberto Pietrini
# Statisticalhitfinder is distributed under the terms of the Simplified BSD License.
# -------------------------------------------------------------------------

import os
os.system("source /home/alberto/.bashrc")

import numpy as np
from utils.functions1 import Functions
from scipy.stats import mode
from utils.parallelH5 import mpi_h5rw
from utils.functions2 import find_1ph_peak, baglivo
from utils.misses import prel_misses
import time
import ConfigParser

# PARSES A CONFIGURATION FILE
""" the configuration file is stored in the folder cfgfile_cxidb-XX/ and contains:
    - experiment name
    - dark run number
    - sample run number (in the case of CXIDB files, only the number as listed in CXI Databank has to be put - 78 in this case - since the CXIDB file is "cxidb-78.cxi")
    - ADU correction to apply (per-ASIC or per-colum per-ASIC)
"""

Config = ConfigParser.ConfigParser()
Config.read("/home/alberto/Statisticalhitfinder/cfgfile/cfgfile_%s.ini" %(os.environ['SLURM_JOB_NAME']))
print "JOB NAME", os.environ['SLURM_JOB_NAME']

def CfgSetMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
        except:
            dict1[option] = None

    return dict1

def dark_mode():
    """ Retrieves the dark mode.

    Reads a dark run file of raw ADU values saved in CXI format, finds the mode in each
    pixel and saves the corresponding ADU values.

    Args:
        None

    Returns:
        Saves to disk the dark mode values of each pixels in HDF5 format. The file created
        is 'dark_mode_back'.
    """

    # READS DARK RUN FILE
    fileID = mpirw.createH5(dir_darkrun + 'run%s_done.h5' %(dark_runnr),'r')
    dark_run = mpirw.read_ds(fileID, 'back')
    dimensions = dark_run.shape

    # WRITES THE FILE WITH THE MODE STORED FOR EACH PIXEL
    mode_fileID = mpirw.createH5(wdir+'dark_mode_back.h5','w')
    dark_mode = mpirw.create_ds(mode_fileID,'mode', (1,dimensions[1],dimensions[2]), 'int16')


    # DARK MODE OF EACH PIXEL
    # Here we use only 5000 frames/events to compute the dark mode
    start, end, ext_dim = dist_load.chunk_file(dimensions,'pixelwise')
    pixels = dark_run[:,start:end,:].reshape([dimensions[0],ext_dim])

    dr_mode = mode(pixels[0:5000,:],0)[0].reshape([end-start,dimensions[2]])

    # SAVES DATA INTO FILE
    dark_mode[:,start:end,:] = dr_mode
    mpirw.ds_flush(mode_fileID)
    comm.barrier()

    # CLOSES FILES
    mpirw.h5_close(fileID)
    mpirw.h5_close(mode_fileID)

def common_mode():
    """ Computes a common mode subtraction to obtain corrected ADU values.

    Reads a dark/background/sample run file of raw ADU values saved in CXI format and
    subtracts the dark common mode. In addition to that, an offset is subtracted in two possible modes:
    - per-ASIC: the ofset of the ASIC is subtracted from the (RAW - DARK_MODE) data;
    - per-column per-ASIC: in addition to the former, also the offset in each pixel column in each ASIC is subtracted.

    Args:
        None

    Returns:
        Saves to disk a 3D array (#events X #pixels_xaxis X #pixels_yaxis) of corrected ADU values in the mode selected, in HDF5
        format. The file created is 'common_mode_back_column' (or 'common_mode_back_offset').
        The mode can be selected in the configuration file.
    """

    # PRINT RANK
    print "RANK - ", rank

    # READS DARK/BACKGROUND/SAMPLE FILE
    sample_fileID = mpirw.createH5(dir_raw_data + 'cxidb-%s.cxi' %(runnr),'r')
    print "FILEID", sample_fileID
    tms = mpirw.read_ds(sample_fileID, 'LCLS/machineTime')
    indxx = np.where(tms[:] != 0)
    sample_run = mpirw.read_ds(sample_fileID, 'entry_1/data_1/data')
    dimensions = (indxx[0].shape[0], sample_run.shape[1], sample_run.shape[2])

    # PRINT AFTER CREATING FILE
    print "File Created - ", rank

    # READS AND WRITES GMD VALUES
    gmd_values = mpirw.read_ds(sample_fileID, 'cheetah/unshared/gmd1') #unshared in cxi version 1.30; event_data in cxi version >=1.40 
    dim_gmd = dimensions[0]
    gmd_fileID = mpirw.createH5(wdir+'gmd.h5','w')
    gmd = mpirw.create_ds(gmd_fileID,'cheetah/event_data/gmd1', (dim_gmd,), 'float32') #event_data
    gmd[:] = gmd_values[:dimensions[0]]

    # READS AND WRITES TIMESTAMPS
    dim_tsp = dimensions[0]
    timestamp_fileID = mpirw.createH5(wdir+'timestamp.h5','w')
    timestamp = mpirw.create_ds(timestamp_fileID,'LCLS/machineTime', (dim_tsp,), 'float64')
    timestamp[:] = tms[:][indxx]

    # WRITES THE COMMON MODE CORRECTED FILE
    common_mode_fileID = mpirw.createH5(wdir+'common_mode_back_%s.h5' %(end_filename),'w')
    common_mode = mpirw.create_ds(common_mode_fileID,'data', dimensions, 'int16')
    common_ramp = mpirw.create_ds(common_mode_fileID,'ramp', dimensions, 'int16')
    
    # GETS ASIC DARK MODE
    #dark_mode_fileID = mpirw.createH5('./files/dark_mode_back.h5','r')
    dark_mode = mpirw.read_ds(sample_fileID, 'entry_1/data_analysis/dark_mode')

    # COMMON MODE 
    start, end = dist_load.chunk_file(dimensions,'framewise')
    frames = sample_run[start:end]

    def split_4(x,a,b):
        """ Perform the ADU correction ASIC by ASIC.

        Args:
            x - dark mode subtracted data
            a - half of the number of pixel along x-axis
            b - half of the number of pixel along y-axis

        Returns:
            corr_frame: corrected ADU values in one of the 2 modes.
            ramp: ADU values of the offsets along the column

            Those returns are saved into the file 'common_mode_back_XXX.h5' [XXX can be 'column' or 'offset']
        """

        #x = x[0] # because 'dark_subtraction' is (1,370,388)
        p1, p2, p3, p4 = x[:a,:b], x[:a,b:], x[a:,b:], x[a:,:b]
        m1p, m2p, m3p, m4p = np.median(p1), np.median(p2), np.median(p3), np.median(p4)
               
        t1, t2, t3, t4 = (p1 - np.tile(m1p,p1.shape),  p2 - np.tile(m2p,p2.shape),
                         p3 - np.tile(m3p,p3.shape), p4 - np.tile(m4p,p4.shape))

        # per-ASIC offset 
        o1, o2, o3, o4 = t1, t2, t3, t4

        hs1 = np.hstack([o1,o2])
        hs2 = np.hstack([o4,o3])
        pA_offset = np.vstack([hs1,hs2])

        # offset along the column (ramp)
        r1, r2, r3, r4 = (np.tile(np.median(t1,0),[t1.shape[0],1]), np.tile(np.median(t2,0),[t2.shape[0],1]),
                         np.tile(np.median(t3,0),[t3.shape[0],1]), np.tile(np.median(t4,0),[t4.shape[0],1]))

        rhs1 = np.hstack([r1,r2])
        rhs2 = np.hstack([r4,r3])
        ramp = np.vstack([rhs1,rhs2])

        # per-column per-ASIC offset
        pCpA_offset = pA_offset - ramp

        if ADUcorr==0:
            corr_frame = pA_offset
        elif ADUcorr==1:
            corr_frame = pCpA_offset
        else:
            print "Incorrect ADUcorr value", ADUcorr

        return corr_frame, ramp

    for i,frame in enumerate(frames):
        dark_subtraction = frame - dark_mode[:] 
        offset_sub, ramps = split_4(dark_subtraction,185,194) #we are here dealiing with cspad2x2 -> (400,400)
        common_mode[start+i] = offset_sub #dark_subtraction - offset - row_mode
        common_ramp[start+i] = ramps
        mpirw.ds_flush(common_mode_fileID)

    comm.barrier()

    # CLOSES FILES
    mpirw.h5_close(sample_fileID)
    mpirw.h5_close(gmd_fileID)
    mpirw.h5_close(timestamp_fileID)
    mpirw.h5_close(common_mode_fileID)
    #mpirw.h5_close(dark_mode_fileID)

def gain():
    # READ COMMON MODE SUBTRACTED FILE
    cms_fileID = mpirw.createH5(wdir+'common_mode_back_%s.h5' %(end_filename),'r')
    cms_run = mpirw.read_ds(cms_fileID, 'data')
    dimensions = cms_run.shape

    # WRITE THE GAINMAP FILE
    gain_fileID = mpirw.createH5(wdir+'gainmap_%s.h5' %(end_filename),'w')
    gain = mpirw.create_ds(gain_fileID,'data', (9, dimensions[1], dimensions[2]), 'float32')

    # FIND THE GAIN FOR EACH PIXEL
    start, end, ext_dim = dist_load.chunk_file(dimensions,'pixelwise')
    pixels = cms_run[:,start:end,:].reshape([dimensions[0],ext_dim])
    single_photon_peak = np.apply_along_axis(find_1ph_peak,0,pixels).reshape([9,end-start,dimensions[2]])

    #(REMOVE)print "Start saving, rank nr.", rank, single_photon_peak.shape

    # SAVE DATA INTO FILE
    gain[:,start:end,:] = single_photon_peak
    mpirw.ds_flush(gain_fileID)
    comm.barrier()

    # CLOSE FILES
    mpirw.h5_close(cms_fileID)
    mpirw.h5_close(gain_fileID)

def photon_count():
    # READ GMD FILE SAMPLE
    gmd_fileID = mpirw.createH5(wdir+'gmd.h5','r')
    gmd_values = mpirw.read_ds(gmd_fileID, 'cheetah/event_data/gmd1')
    gmd_sample = gmd_values[:]
    gmd_bkg = gmd_sample

    # READ COMMON MODE SUBTRACTED FILE
    cms_fileID = mpirw.createH5(wdir+'common_mode_back_%s.h5' %(end_filename),'r')
    cms_run = mpirw.read_ds(cms_fileID, 'data')
    dimensions = cms_run.shape

    # READ THE GAINMAP FILE
    gain_fileID = mpirw.createH5(dir_raw_data + 'cxidb-%s.cxi' %(runnr),'r')
    gain = mpirw.read_ds(gain_fileID, 'entry_1/data_analysis/gainmap')
    gainmap = gain


    # WRITE THE PHOTON COUNT FILE
    photon_space_frames_fileID = mpirw.createH5(wdir+'photon_space_frames_%s.h5' %(end_filename),'w')
    ph_space_ds = mpirw.create_ds(photon_space_frames_fileID,'data', dimensions, 'int16')
    ph_count_per_frame = mpirw.create_ds(photon_space_frames_fileID,'photon_count_per_frame', (dimensions[0],), 'int32')
    ph_sum_over_frames_worker = mpirw.create_ds(photon_space_frames_fileID,'sum_over_frames_worker', (nr_workers,dimensions[1],dimensions[2]), 'int64')

    # index of frames to be treated as background
    index_bkg = mpirw.create_ds(photon_space_frames_fileID,'indices1', (dimensions[0],), 'int8')

    # index where you exclude the extreme outliers - above and below 2 sigma
    index_up_down = mpirw.create_ds(photon_space_frames_fileID,'indices2', (dimensions[0],), 'int8')

    # fitted gmd values - background and sample
    fitv_bkg = mpirw.create_ds(photon_space_frames_fileID,'fitv_bkg', (dimensions[0],), 'float32')
    fitv_sample = mpirw.create_ds(photon_space_frames_fileID,'fitv_sample', (dimensions[0],), 'float32')

    # LOAD POISSON MASK
    try:
        mask_fileID = mpirw.createH5(wdir+'poisson_mask.h5','r')
        mask = mpirw.read_ds(mask_fileID, 'data')[:].astype('bool')
        print "Mask found!"
    except Exception:
        mask = np.ones([dimensions[1],dimensions[2]]).astype('bool')

    # SAVE PHOTON COUNT FRAMES INTO FILES
    start, end = dist_load.chunk_file(dimensions,'framewise')
    frames = cms_run[start:end]

    for i,frame in enumerate(frames):
        ph_space_fr = (frame/gainmap).round()
        ph_space_fr[ph_space_fr<0] = 0
        ph_space_fr[np.isnan(ph_space_fr)] = 0

        ph_space_ds[start+i] = ph_space_fr
        ph_count_per_frame[start+i] = ph_space_fr[mask].sum()
    mpirw.ds_flush(photon_space_frames_fileID)
    comm.barrier()

    prel_ms = prel_misses(ph_count_per_frame[start:end],gmd_bkg[start:end],mask)
    comm.barrier()
    fake_bkg = prel_ms

    index_bkg[start:end] = fake_bkg.astype(int)
    index_up_down[start:end] = fake_bkg.astype(int)
    ph_sum_over_frames_worker[rank] = np.sum(ph_space_ds[start:end][fake_bkg],0)
    mpirw.ds_flush(photon_space_frames_fileID)
    comm.barrier()

    phcn = ph_count_per_frame[:][index_up_down[:].astype(bool)]
    gmdn = gmd_bkg[index_up_down[:].astype(bool)]

    p3 = np.polyfit(gmdn,phcn,3)
    fit_bkg = np.dot(p3,np.array([gmd_bkg**3,gmd_bkg**2,gmd_bkg,1]))
    fit_sample = np.dot(p3,np.array([gmd_sample**3,gmd_sample**2,gmd_sample,1]))

    ## SAVE THE FITTED VALUES   
    fitv_bkg[:] = fit_bkg
    fitv_sample[:] = fit_sample
    mpirw.ds_flush(photon_space_frames_fileID)

    comm.barrier()
    print "File saved - rank nr.", rank

    # CLOSE FILES
    mpirw.h5_close(gmd_fileID)
    mpirw.h5_close(cms_fileID)
    mpirw.h5_close(gain_fileID)
    mpirw.h5_close(photon_space_frames_fileID)
    try:
        mpirw.h5_close(mask_fileID)
    except Exception:
        pass

def lambda_values():
    # READ PHOTON COUNT FILE OF SAMPLE
    phot_fileID = mpirw.createH5(wdir+'photon_space_frames_%s.h5' %(end_filename),'r')
    phot_run = mpirw.read_ds(phot_fileID, 'data')
    phot_per_frame = mpirw.read_ds(phot_fileID, 'photon_count_per_frame')
    phot_per_frame_sample = phot_per_frame[:]
    phot_per_frame_bkg = phot_per_frame[:]
    dimensions = phot_run.shape

    # READ PHOTON COUNT FILE OF SAMPLE
    sum_over_bkg = mpirw.read_ds(phot_fileID, 'sum_over_frames_worker')
    matrix_sum = np.sum(sum_over_bkg[:],0)

    # WRITE THE LAMBDA MATRIX
    lambdas_fileID = mpirw.createH5(wdir+'lambda_matrix_column.h5','w')
    lambdas_ds = mpirw.create_ds(lambdas_fileID,'data', dimensions, 'float32')

    # SAVE PHOTON COUNT FRAMES INTO FILES
    start, end = dist_load.chunk_file(dimensions,'framewise')

    # FITTING GMD VALUES TO PHOTON COUNT
    fit_sample = mpirw.read_ds(phot_fileID, 'fitv_sample')[:]
    index_bkg = mpirw.read_ds(phot_fileID,'indices1')[:]
    fit_bkg = mpirw.read_ds(phot_fileID, 'fitv_bkg')[:][index_bkg.astype(bool)]

    for i,el in enumerate(fit_sample[start:end]):
        lambda_frame = el*matrix_sum/fit_bkg.sum()
        lambda_frame[lambda_frame<=0] = 1e-30
        lambdas_ds[start+i] = lambda_frame
        mpirw.ds_flush(lambdas_fileID)
    comm.barrier()
    print "File saved - rank nr.", rank

    # CLOSE FILES
    mpirw.h5_close(phot_fileID)
    mpirw.h5_close(lambdas_fileID)

def poissmask():
    # READ PHOTON COUNT FILE OF SAMPLE
    phot_fileID = mpirw.createH5(wdir+'photon_space_frames_%s.h5' %(end_filename),'r')
    
    print "ID file", phot_fileID
    phot_run = mpirw.read_ds(phot_fileID, 'data')
    phot_per_frame = mpirw.read_ds(phot_fileID, 'photon_count_per_frame')
    phot_per_frame_sample = phot_per_frame[:]
    phot_per_frame_bkg = phot_per_frame[:]
    dimensions = (40000, phot_run.shape[1], phot_run.shape[2])
    print "shape", dimensions

    # READ LAMBDAS
    lambdas_fileID = mpirw.createH5(wdir+'lambda_matrix_%s.h5' %(end_filename),'r')
    lambdas_ds = mpirw.read_ds(lambdas_fileID, 'data')

    print "Lambdas read"

    # WRITE THE POISSON MASK
    pmsk_fileID = mpirw.createH5(wdir+'poisson_mask.h5','w')
    print "mask_ID", pmsk_fileID
    pmsk_ds = mpirw.create_ds(pmsk_fileID,'data', (dimensions[1],dimensions[2]), 'int16')
    print "data2mask", pmsk_ds, pmsk_ds.shape

    # SAVE PHOTON COUNT FRAMES INTO FILES
    start, end, ext_dim = dist_load.chunk_file(dimensions,'pixelwise')

    print "Chunks'dimension", start, end, ext_dim
    # FITTING GMD VALUES TO PHOTON COUNT
    index_bkg = mpirw.read_ds(phot_fileID,'indices1')[:]
    #index = np.where(index_bkg==1)
    index = index_bkg==1
    index = index[0:40000]
    index_shape = np.where(index==True)[0].shape[0]
    
    # SAMPLE LAMBDAS AND COMPARE WITH PHOTONS
    photons = phot_run[index,start:end,:].reshape([index_shape,ext_dim]).T

    print "Photon done", rank
    mpirw.ds_flush(phot_fileID)
    comm.barrier()

    lambdas = lambdas_ds[index,start:end,:].reshape([index_shape,ext_dim]).T
    print "Lambdas", rank
    mpirw.ds_flush(lambdas_fileID)
    comm.barrier()

    sampled_l = np.random.poisson(lambdas)

    values = []
    for i in range(0,ext_dim):
        y1,x1 = np.histogram(photons[i],bins=np.arange(photons[i].min(),photons[i].max()+1))
        y2,x2 = np.histogram(sampled_l[i],bins=np.arange(photons[i].min(),photons[i].max()+1))
        value = np.dot(y1,y2)/(np.sqrt(np.dot(y1,y1))*np.sqrt(np.dot(y2,y2)))
        values.append(value)

    thrs = 0.9999

    vals = np.array(values)
    vals[np.isnan(vals)] = 0
    vals[vals<thrs] = 0
    vals[vals>=thrs] = 1

    vals_res = vals.reshape([end-start, dimensions[2]])
    pmsk_ds[start:end] = vals_res
    mpirw.ds_flush(pmsk_fileID)
    comm.barrier()

    # CLOSE FILES
    mpirw.h5_close(phot_fileID)
    mpirw.h5_close(lambdas_fileID)
    mpirw.h5_close(pmsk_fileID)


def baglivo_score():
    # READ PHOTON COUNT FILE FROM BACKGROUND
    photon_space_frames_fileID = mpirw.createH5(wdir+'photon_space_frames_%s.h5' %(end_filename),'r')
    ph_space_ds = mpirw.read_ds(photon_space_frames_fileID, 'data')
    fit_sample = mpirw.read_ds(photon_space_frames_fileID, 'fitv_sample')
    dimensions = ph_space_ds.shape

    # READ LAMBDAS
    lambdas_fileID = mpirw.createH5(wdir+'lambda_matrix_%s.h5' %(end_filename),'r')
    lambdas_ds = mpirw.read_ds(lambdas_fileID, 'data')

    # WRITE BAGLIVO SCORE FILE
    score_fileID = mpirw.createH5(wdir+'baglivo_score.h5','w')
    score_ds = mpirw.create_ds(score_fileID,'data', (dimensions[0],), 'float32')

    # LOAD MASK
    try:
        mask_fileID = mpirw.createH5(wdir+'poisson_mask.h5','r')
        mask = mpirw.read_ds(mask_fileID, 'data')[:].astype('bool')
        print "Mask found!"
    except Exception:
        mask = np.ones([dimensions[1],dimensions[2]]).astype('bool')


    # SAVE SCORES FRAMES INTO FILES
    start, end = dist_load.chunk_file(dimensions,'framewise')
    for i,el in enumerate(fit_sample[start:end]):
        score_ds[start+i] = baglivo(lambdas_ds[start+i][mask].ravel(),ph_space_ds[start+i][mask].ravel(),el)

        if el<=400: score_ds[start+i] = 0

        mpirw.ds_flush(score_fileID)

    comm.barrier()

    # CLOSE FILES
    mpirw.h5_close(photon_space_frames_fileID)
    mpirw.h5_close(lambdas_fileID)
    mpirw.h5_close(score_fileID)
    try:
        mpirw.h5_close(mask_fileID)
    except Exception:
        pass


def main():
    """If one of the steps has already been made and there is no need to re-run it, just comment it out"""
    #dark_mode()
    #common_mode()
    #gain()
    photon_count()
    lambda_values()
    poissmask()
    photon_count()
    lambda_values()
    baglivo_score()

if __name__ == '__main__':
    mpirw = mpi_h5rw()
    nr_workers, rank, comm = mpirw.mpi_params()
    dist_load = Functions(nr_workers,rank)

    nameExp = CfgSetMap("exp_details")['exp']
    dark_runnr = np.int16(CfgSetMap("exp_details")['dark_runnr'])
    runnr = np.int16(CfgSetMap("exp_details")['sample_runnr'])

    # BACK DETECTOR
    #dir_darkrun = '/scratch/fhgfs/LCLS/cxi/%s/hdf5/' %(nameExp)
    dir_raw_data = '/scratch/fhgfs/alberto/'
    wdir = '/scratch/fhgfs/alberto/testPROVA/'

    if rank==0:
        if not os.path.isdir(wdir):
            os.mkdir(wdir)
        else:
            pass

    comm.barrier()

    # define ADUcorr
    ADUcorr = np.int16(CfgSetMap("common_mode")['adu_corr'])

    if ADUcorr == 0:
        end_filename = 'offset'
    elif ADUcorr == 1:
        end_filename = 'column'

    if rank == 0: t0 = time.time()

    main()
    mpirw.comm.barrier()

    if rank == 0:
        tf = time.time()
        print "Writing done! Elapsed time:", tf - t0, "seconds."

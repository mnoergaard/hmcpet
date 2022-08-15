#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: martinnorgaard
"""

import os
import sys
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl as fsl
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataSink
from nipype import Node, Function, MapNode
from nipype.pipeline import Workflow
from nipype.interfaces.io import SelectFiles
from bids import BIDSLayout

# Define BIDS directory, including input/output relations
def main():
    """
    Main function to run the workflow.

    Arguments
    ---------
    bids_dir : BIDS directory
    n_procs : number of processors to use
    """

    args = sys.argv[1:]
    bids_dir = args[0]
    n_procs = int(args[1])

    layout = BIDSLayout(bids_dir)
    infosource = Node(IdentityInterface(
                    fields = ['subject_id','session_id']),
                    name = "infosource")

    infosource.iterables = [('subject_id', layout.get_subjects()), 
                        ('session_id', layout.get_sessions())]

    templates = {'pet': 'sub-{subject_id}/ses-{session_id}/pet/*_pet.[n]*', 
             'json': 'sub-{subject_id}/ses-{session_id}/pet/*_pet.json'}
           
    selectfiles = Node(SelectFiles(templates, 
                               base_directory = bids_dir), 
                               name = "select_files")

    datasink = Node(DataSink(base_directory = os.path.join(bids_dir, 'derivatives')), 
                         name = "datasink")

    substitutions = [('_subject_id', 'sub'), ('_session_id_', 'ses')]
    subjFolders = [('sub-%s_ses-%s' % (sub,ses), 'sub-%s/ses-%s' %(sub,ses))
               for ses in layout.get_sessions()
               for sub in layout.get_subjects()]

    substitutions.extend(subjFolders)
    datasink.inputs.substitutions = substitutions

    # Define nodes for hmc workflow

    split_pet = Node(interface = fs.MRIConvert(split = True), 
                     name = "split_pet")
    
    smooth_frame = MapNode(interface=fsl.Smooth(fwhm=10), 
                           name="smooth_frame", 
                           iterfield=['in_file'])
    
    thres_frame = MapNode(interface = fsl.maths.Threshold(thresh = 20, use_robust_range = True),
                          name = "thres_frame", 
                          iterfield = ['in_file'])
    
    estimate_motion = Node(interface = fs.RobustTemplate(auto_detect_sensitivity = True,
                                            intensity_scaling = True,
                                            average_metric = 'mean',
                                            args = '--cras'),
                           name="estimate_motion", iterfield=['in_files'])
    
    correct_motion = MapNode(interface = fs.ApplyVolTransform(), 
                             name = "correct_motion", 
                             iterfield = ['source_file', 'reg_file', 'transformed_file'])
    
    concat_frames = Node(interface = fs.Concatenate(concatenated_file = 'mc.nii.gz'), 
                         name = "concat_frames")
    
    lta2xform = MapNode(interface = fs.utils.LTAConvert(), 
                        name = "lta2xform", 
                        iterfield = ['in_lta', 'out_fsl'])
    
    est_trans_rot = MapNode(interface = fsl.AvScale(all_param = True), 
                            name = "est_trans_rot", 
                            iterfield = ['mat_file'])
    
    est_min_frame = Node(Function(input_names = ['json_file'],
                                  output_names = ['min_frame'],
                                  function = get_min_frame),
                         name = "est_min_frame")
    
    upd_frame_list = Node(Function(input_names = ['in_file','min_frame'],
                                   output_names = ['upd_list_frames'],
                                   function = update_list_frames),
                          name = "upd_frame_list")
    
    upd_transform_list = Node(Function(input_names = ['in_file','min_frame'],
                                       output_names = ['upd_list_transforms'],
                                       function = update_list_transforms),
                          name = "upd_transform_list")
    
    hmc_movement_output = Node(Function(input_names = ['translations', 'rot_angles', 'rotation_translation_matrix','in_file'],
                                           output_names = ['hmc_confounds'],
                                           function = combine_hmc_outputs),
                               name = "hmc_movement_output")
    
    plot_motion = Node(Function(input_names = ['in_file'],
                                           function = plot_motion_outputs),
                               name = "plot_motion")
    
# Connect workflow - init_pet_hmc_wf
    workflow = Workflow(name = "hmc_workflow")
    workflow.config['execution']['remove_unnecessary_outputs'] = 'false'
    workflow.base_dir = bids_dir
    workflow.connect([(infosource, selectfiles, [('subject_id', 'subject_id'),('session_id', 'session_id')]), 
                         (selectfiles, split_pet, [('pet', 'in_file')]),
                         (selectfiles, est_min_frame, [('json', 'json_file')]),
                         (split_pet,smooth_frame,[('out_file', 'in_file')]),
                         (smooth_frame,thres_frame,[('smoothed_file', 'in_file')]),
                         (thres_frame,upd_frame_list,[('out_file', 'in_file')]),
                         (est_min_frame,upd_frame_list,[('min_frame', 'min_frame')]),
                         (upd_frame_list,estimate_motion,[('upd_list_frames', 'in_files')]),
                         (thres_frame,upd_transform_list,[('out_file', 'in_file')]),
                         (est_min_frame,upd_transform_list,[('min_frame', 'min_frame')]),
                         (upd_transform_list,estimate_motion,[('upd_list_transforms', 'transform_outputs')]),
                         (split_pet,correct_motion,[('out_file', 'source_file')]),
                         (estimate_motion,correct_motion,[('transform_outputs', 'reg_file')]),
                         (estimate_motion,correct_motion,[('out_file', 'target_file')]),
                         (split_pet,correct_motion,[(('out_file', add_mc_ext), 'transformed_file')]),
                         (correct_motion,concat_frames,[('transformed_file', 'in_files')]),
                         (estimate_motion,lta2xform,[('transform_outputs', 'in_lta')]),
                         (estimate_motion,lta2xform,[(('transform_outputs', lta2mat), 'out_fsl')]),
                         (lta2xform,est_trans_rot,[('out_fsl', 'mat_file')]),
                         (est_trans_rot,hmc_movement_output,[('translations', 'translations'),('rot_angles', 'rot_angles'),('rotation_translation_matrix','rotation_translation_matrix')]),
                         (upd_frame_list,hmc_movement_output,[('upd_list_frames', 'in_file')]),
                         (hmc_movement_output,plot_motion,[('hmc_confounds','in_file')])
                         ])
    wf = workflow.run(plugin='MultiProc', plugin_args={'n_procs' : n_procs})

# HELPER FUNCTIONS
def update_list_frames(in_file, min_frame):   
    """  
    Function to update the list of frames to be used in the hmc workflow.

    Parameters
    ----------
    in_file : list of frames
    min_frame : minimum frame to use for the analysis (first frame after 2 min)

    Returns
    -------
    new_list : list of updated frames to be used in the hmc workflow

    """
    
    new_list = [in_file[min_frame]] * min_frame + in_file[min_frame:]
    return new_list

def update_list_transforms(in_file, min_frame):   
    """
    Function to update the list of transforms to be used in the hmc workflow.
    
    Parameters
    ----------
    in_file : list of transforms
    min_frame : minimum frame to use for the analysis (first frame after 2 min)

    Returns
    -------
    lta_list : list of updated transforms to be used in the hmc workflow
    """
    
    new_list = [in_file[min_frame]] * min_frame + in_file[min_frame:]
    lta_list = [ext.replace('nii.gz','lta') for ext in new_list]  
    return lta_list

def add_mc_ext(in_file):    
    """ 
    Function to add the mc extension to the list of file names.

    Parameters
    ----------
    in_file : file name to be updated

    Returns
    -------
    mc_list : list of updated file names with mc extension
    """
    
    if len(in_file) > 1:
        mc_list = [ext.replace('.nii.gz','_mc.nii') for ext in in_file] # list comphrehension
    else:
        mc_list = in_file.replace('.nii.gz','_mc.nii')
    return mc_list

def lta2mat(in_file):  
    """
    Function to convert the lta file to the fsl format (.mat).
    
    Parameters
    ----------
    in_file : list of lta files to be converted

    Returns
    -------
    mat_list : list of mat files
    """
    
    mat_list = [ext.replace('.lta','.mat') for ext in in_file]
    return mat_list 

def get_min_frame(json_file):  
    """
    Function to extract the frame number after 120 seconds of mid frames of dynamic PET data to be used with motion correction
        
    Parameters
    ----------
    json_file : json file containing the frame length and duration of the dynamic PET data

    Returns
    -------
    min_frame : minimum frame to use for the motion correction (first frame after 2 min)
    """  
    
    import os
    from os.path import join, isfile
    import numpy as np
    import json

    with open(json_file, 'r') as f:
        info = json.load(f)
        frames_duration = np.array(info['FrameDuration'], dtype=float)
        frames_start =np.pad(np.cumsum(frames_duration)[:-1],(1,0))
        mid_frames = frames_start + frames_duration/2

        min_frame = next(x for x, val in enumerate(mid_frames)
                                  if val > 120)    
    return min_frame  


def combine_hmc_outputs(translations, rot_angles, rotation_translation_matrix, in_file):   
    """
    
    Function to combine the outputs of the hmc workflow.

    Parameters
    ----------
    translations : list of estimated translations across frames
    rot_angles : list of estimated rotational angles across frames
    rotation_translation_matrix : list of estimated rotation translation matrices across frames
    in_file : list of frames to be used in the hmc workflow

    Returns
    -------
    Output path to confounds file for head motion correction
    """
    
    import os
    import pandas as pd
    import numpy as np
    import nibabel as nib
    
    new_pth = os.getcwd()
    
    movement = []
    for idx, trans in enumerate(translations):
        
        img = nib.load(in_file[idx])
        vox_ind = np.asarray(np.nonzero(img.get_fdata()))
        pos_bef = np.concatenate((vox_ind,np.ones((1,len(vox_ind[0,:])))))
        pos_aft = rotation_translation_matrix[idx] @ pos_bef
        diff_pos = pos_bef-pos_aft
        
        max_x = abs(max(diff_pos[0,:], key=abs))
        max_y = abs(max(diff_pos[1,:], key=abs))
        max_z = abs(max(diff_pos[2,:], key=abs))
        overall = np.sqrt(diff_pos[0,:] ** 2 + diff_pos[1,:] ** 2 + diff_pos[2,:] ** 2)
        max_tot = np.max(overall)
        median_tot = np.median(overall)
        
        movement.append(np.concatenate((translations[idx],rot_angles[idx], [max_x, max_y, max_z, max_tot, median_tot])))
        
    confounds = pd.DataFrame(movement, columns=['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'max_x', 'max_y', 
                                                'max_z', 'max_tot', 'median_tot'])
    confounds.to_csv(os.path.join(new_pth,'hmc_confounds.tsv'), sep='\t')
    # np.savetxt(os.path.join(new_pth,'hmc_confounds.tsv'), movement, fmt='%10.5f', delimiter='\t', header='trans_x    trans_y    trans_z    rot_x    rot_y    rot_z')
    
    return os.path.join(new_pth,'hmc_confounds.tsv')

def plot_motion_outputs(in_file):   
    """
    Function to plot estimated motion data
    
    Parameters 
    ----------
    in_file : list of estimated motion data

    Returns
    -------
    Plots of estimated motion data
    """
    
    import os
    import pandas as pd
    import numpy as np
    import nibabel as nib
    import matplotlib.pyplot as plt
    
    confounds = pd.read_csv(in_file, sep='\t')
    
    new_pth = os.getcwd()
    
    n_frames = len(confounds.index)
    
    plt.figure(figsize=(11,5))
    plt.plot(np.arange(0,n_frames), confounds['trans_x'], "-r", label='trans_x')
    plt.plot(np.arange(0,n_frames), confounds['trans_y'], "-g", label='trans_y')
    plt.plot(np.arange(0,n_frames), confounds['trans_z'], "-b", label='trans_z')
    plt.legend(loc="upper left")
    plt.ylabel('Translation [mm]')
    plt.xlabel('frame #')
    plt.grid(visible=True)
    plt.savefig(os.path.join(new_pth,'translation.png'), format='png')
    plt.close()
    
    plt.figure(figsize=(11,5))
    plt.plot(np.arange(0,n_frames), confounds['rot_x'], "-r", label='rot_x')
    plt.plot(np.arange(0,n_frames), confounds['rot_y'], "-g", label='rot_y')
    plt.plot(np.arange(0,n_frames), confounds['rot_z'], "-b", label='rot_z')
    plt.legend(loc="upper left")
    plt.ylabel('Rotation [degrees]')
    plt.xlabel('frame #')
    plt.grid(visible=True)
    plt.savefig(os.path.join(new_pth,'rotation.png'), format='png')
    plt.close()
    
    plt.figure(figsize=(11,5))
    plt.plot(np.arange(0,n_frames), confounds['max_x'], "--r", label='max_x')
    plt.plot(np.arange(0,n_frames), confounds['max_y'], "--g", label='max_y')
    plt.plot(np.arange(0,n_frames), confounds['max_z'], "--b", label='max_z')
    plt.plot(np.arange(0,n_frames), confounds['max_tot'], "-k", label='max_total')
    plt.plot(np.arange(0,n_frames), confounds['median_tot'], "-m", label='median_tot')
    plt.legend(loc="upper left")
    plt.ylabel('Movement [mm]')
    plt.xlabel('frame #')
    plt.grid(visible=True)
    plt.savefig(os.path.join(new_pth,'movement.png'), format='png')
    plt.close()
    
if __name__ == '__main__':
    main()
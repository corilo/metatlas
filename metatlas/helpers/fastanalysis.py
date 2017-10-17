import sys, os

from metatlas.helpers import metatlas_get_data_helper_fun as ma_data
from metatlas.helpers import dill2plots as dp
from metatlas.helpers import chromplotplus as cpp
from metatlas.helpers import spectralprocessing as sp

import numpy as np
import multiprocessing as mp
import pandas as pd
import copy


def scores_for_each_compound(atlas_df, metatlas_dataset):
    """
    Returns pandas dataframe with columns 'max_intensity', 'median_rt_shift','median_mz_ppm', 'max_msms_score', 
    'num_frag_matches', and 'max_relative_frag_intensity', rows of compounds in metatlas_dataset, and values
    of the best "score" for a given compound across all files. 
    
    'max_intensity': highest EIC across all files for given compound
    'median_rt_shift': shift of median RT across all files for given compound to reference
    'median_mz_ppm': ppm of median mz across all files for given compound relative to reference
    'max_msms_score': highest compound dot-product score across all files for given compound relative to reference
    'num_frag_matches': number of matching mzs when calculating max_msms_score
    'max_relative_frag_intensity': ratio of second highest to first highest intensity of matching sample mzs


    :param metatlas_dataset:
    :param atlas_df:

    :return scores_df: pandas dataframe
    """
    
    file_names = ma_data.get_file_names(metatlas_dataset)
    compound_names = ma_data.get_compound_names(metatlas_dataset)[0]
    frag_refs = pd.read_json(os.path.join('/project/projectdirs/metatlas/projects/sharepoint/frag_refs.json'))
    
    scores_df = pd.DataFrame(np.nan,
                             columns = ['max_intensity',
                                        'median_rt_shift',
                                        'median_mz_ppm',
                                        'max_msms_score',
                                        'num_frag_matches',
                                        'max_relative_frag_intensity'],
                             index = atlas_df.label)
    
    for compound_idx in range(len(compound_names)):
        #max intensity
        max_intensity = np.nan
        for file_idx in range(len(file_names)):
            if len(metatlas_dataset[file_idx][compound_idx]['data']['eic']['intensity']) == 0:
                continue
            file_max_intensity = max(metatlas_dataset[file_idx][compound_idx]['data']['eic']['intensity'])
            if file_max_intensity > max_intensity or np.isnan(max_intensity):
                max_intensity = file_max_intensity
        
        #median rt shift
        compound_ref_rt_peak = metatlas_dataset[0][compound_idx]['identification'].rt_references[0].rt_peak

        median_rt_shift = np.nanmedian([abs(compound_ref_rt_peak - metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['rt_centroid'])
                                       for file_idx in range(len(file_names))])# if isinstance(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['rt_centroid'], float)])
        
        scores_df.iloc[compound_idx].median_rt_shift = median_rt_shift
        
        #median mz ppm
        compound_ref_mz = metatlas_dataset[0][compound_idx]['identification'].mz_references[0].mz

        median_mz_ppm = np.nanmedian([1e6*(abs(compound_ref_mz - metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['mz_centroid']) / compound_ref_mz)
                                      for file_idx in range(len(file_names))])# if isinstance(metatlas_dataset[file_idx][compound_idx]['data']['ms1_summary']['mz_centroid'], float)])
        
        #max msms score
        file_idx, max_msms_score, msv_ref = dp.file_with_max_score(metatlas_dataset, frag_refs, compound_idx, 'inchi_key and polarity')
        
        #number of frag matches and maximum relative frag intensity
        num_frag_matches = np.nan
        max_relative_frag_intensity = np.nan
        if ~np.isnan(max_msms_score):
            msv_ref = np.array(msv_ref[0]).T
            msv_sample = sp.sort_ms_vector_by_mz(np.array([metatlas_dataset[file_idx][compound_idx]['data']['msms']['data']['mz'], metatlas_dataset[file_idx][compound_idx]['data']['msms']['data']['i']]))
            
            msv_sample_matches = sp.partition_ms_vectors(msv_sample, msv_ref, .005, 'shape')[0]
            num_frag_matches = len(msv_sample_matches[0])
            
            if num_frag_matches > 1:
                msv_sample_matches_by_intensity = msv_sample_matches[:, msv_sample_matches[1].argsort()]

                max_relative_frag_intensity = msv_sample_matches_by_intensity[1,-2] / msv_sample_matches_by_intensity[1,-1]
            
            
        #assign scores
        scores_df.iloc[compound_idx] = (max_intensity, median_rt_shift, median_mz_ppm, max_msms_score, num_frag_matches, max_relative_frag_intensity)

    return scores_df
    
def filter_metatlas_dataset_by_scores(scores_df, metatlas_dataset, min_intensity, rt_tolerance, mz_tolerance, min_msms_score, allow_no_msms, min_num_frag_matches, min_relative_frag_intensity, full_remove = False):
    
    file_names = ma_data.get_file_names(metatlas_dataset)
    
    compounds_to_remove = []
    
    for compound_idx,(label,row) in enumerate(scores_df.iterrows()):
        if row.max_intensity > min_intensity:
            if row.median_rt_shift < rt_tolerance:
                if row.median_mz_ppm < mz_tolerance:
                    if (row.max_msms_score > min_msms_score) or (allow_no_msms and np.isnan(row.max_msms_score)):
                        if row.num_frag_matches > min_num_frag_matches:
                            if row.max_relative_frag_intensity > min_relative_frag_intensity:
                                continue
        compounds_to_remove.append(compound_idx)
        
    for compound_idx in compounds_to_remove[::-1]:
        for file_idx in range(len(file_names)):
            if full_remove:
                metatlas_dataset[file_idx].pop(compound_idx)
            else:
                metatlas_dataset[file_idx][compound_idx] = empty_data

    return metatlas_dataset
    

def filter_and_dump(atlas, groups, output_dir, 
                    min_intensity = 1e5, 
                    rt_tolerance = .25, 
                    mz_tolerance = 10, 
                    min_msms_score = .7, allow_no_msms = False,
                    min_num_frag_matches = 3,  min_relative_frag_intensity = .05,
                    compress = False):
    
    """
    Creates error bars, chromatograms, and identification figures in output_dir for compounds in
    metatlas_dataset created by atlas and groups which meet the minimum requirements set by
    'min_intensity', 'rt_tolerance','mz_tolerance', 'min_msms_score', 
    'min_num_frag_matches', and 'min_relative_frag_intensity'.
    
    'min_intensity' =< highest EIC across all files for given compound
    'rt_tolerance' >= shift of median RT across all files for given compound to reference
    'mz_tolerance' >= ppm of median mz across all files for given compound relative to reference
    'min_msms_score' =< highest compound dot-product score across all files for given compound relative to reference
    'min_num_frag_matches' =< number of matching mzs when calculating max_msms_score
    'min_relative_frag_intensity' =< ratio of second highest to first highest intensity of matching sample mzs


    :param atlas:
    :param groups:
    :param output_dir:
    """
    
    atlas_df = ma_data.make_atlas_df(atlas)
    atlas_df['label'] = [cid.name for cid in atlas.compound_identifications]
    
    #make metatlas dataset   
    all_files = []
    for my_group in groups:
        for my_file in my_group.items:
            all_files.append((my_file , my_group, atlas_df, atlas))
    pool = mp.Pool(processes=min(8, len(all_files)))
    metatlas_dataset = pool.map(ma_data.get_data_for_atlas_df_and_file, all_files)
    pool.close()
    pool.terminate()
        
    #scores compounds in metatlas dataset
    scores_df = scores_for_each_compound(atlas_df, metatlas_dataset)
    
    #filter dataset by scores
    filtered_dataset = copy.deepcopy(metatlas_dataset)
    filtered_dataset = filter_metatlas_dataset_by_scores(scores_df, filtered_dataset, 
                                                         min_intensity, 
                                                         rt_tolerance, 
                                                         mz_tolerance, 
                                                         min_msms_score, allow_no_msms,
                                                         min_num_frag_matches,
                                                         min_relative_frag_intensity,
                                                         full_remove=True)
    
    #Scores dataframe
    scores_df.to_csv(os.path.join(output_dir,'compound_scores.csv'), sep='\t')
    
    #Chromatograms
    group = 'sort' # 'page' or 'index' or 'sort' or None
    save = True
    share_y = True

    file_names = ma_data.get_file_names(filtered_dataset)
    compound_names = ma_data.get_compound_names(filtered_dataset)[0]
    args_list = []

    chromatogram_str = 'compound_chromatograms'

    if not os.path.exists(os.path.join(output_dir,chromatogram_str)):
        os.makedirs(os.path.join(output_dir,chromatogram_str))

    for compound_idx, my_compound in enumerate(compound_names):
        my_data = list()
        for file_idx, my_file in enumerate(file_names):
            my_data.append(filtered_dataset[file_idx][compound_idx])
        kwargs = {'data': my_data,
                 'file_name': os.path.join(output_dir, chromatogram_str, my_compound+'.pdf'),
                 'group': group,
                 'save': save,
                 'share_y': share_y,
                 'names': file_names}
        args_list.append(kwargs)

    max_processes = 8
    pool = mp.Pool(processes=min(max_processes, len(filtered_dataset[0])))
    pool.map(cpp.chromplotplus, args_list)
    pool.close()
    pool.terminate()
    
    #Error bars
    peak_height = dp.make_output_dataframe(input_fname = '',input_dataset = filtered_dataset, include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_height')
    dp.plot_errorbar_plots(peak_height, output_loc=os.path.join(output_dir,'error_bar_peak_height'))
    
    #Identification figures
    dp.make_identification_figure_v2(input_dataset = filtered_dataset, input_fname = my_file, include_lcmsruns = [],exclude_lcmsruns = [], output_loc=os.path.join(output_dir,'identification'))
   
#    if compress:
#    #Compress
#        timestr = time.strftime("%Y%m%d-%H%M%S")
#        tarball_name = timestr + '_' + os.path.basename(os.path.normpath(output_dir)) + '.tar.gz'
#        %system tar -zcf $tarball_name -C $output_dir .
#    print 'done'
    
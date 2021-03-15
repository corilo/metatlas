import sys, os
from metatlas.tools import fastanalysis as fa
from metatlas.plots import dill2plots as dp
from metatlas.io import metatlas_get_data_helper_fun as ma_data
from metatlas.plots import chromatograms_mp_plots as cp
from metatlas.plots import chromplotplus as cpp
from metatlas.datastructures import metatlas_objects as metob
import qgrid
from ipywidgets import interact, interactive, fixed
import ipywidgets as widgets
from IPython.display import display
import time
import pickle
import dill
import multiprocessing as mp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plot ## from RT -adjustor notebook: 'import matplotlib.pyplot as plt'
import operator
from importlib import reload

from IPython.core.display import Markdown, display, clear_output, HTML
      

# DISPLAY
from  IPython.core.display  import  display, HTML 
display(HTML("<style>.container { width:100% !important; }</style>"))
pd.set_option('display.max_rows', 5000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', 100)




# WRAPPER FUNCTIONS
def setup_notebook(project_directory,project_name,analysis_round):
    """
    sets up output directories

    """
    # DIRECTORIES
    output_subfolder=os.path.join(project_name,analysis_round) 
    output_dir = os.path.join(project_directory,output_subfolder)

    if not os.path.exists(project_directory):
        os.makedirs(project_directory)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("all output files will be written to: "+output_dir)
        
    return(output_dir)

def printmd(string):
    display(Markdown(string))


def setup_rt_adjustment(project_name,my_id,output_dir,polarity="POS",QC_template_filename ='HILICz150_ANT20190824_TPL_QCv3_Unlab_POS'):
    
    ### FIND FILES
    files = dp.get_metatlas_files(experiment = project_name,name = '%AK%',most_recent = True)
    df = metob.to_dataframe(files)
    my_grid = qgrid.QGridWidget(df=df[['experiment','name','username','acquisition_time']])
    print(len(files)+" files found for this experiment")
    
    controlled_vocab = ['QC','InjBl','ISTD'] #add _ to beginning. It will be stripped if at begining  SIG: should we add InjBL? what about extraction control and blanks?
    version_identifier = my_id
    file_dict = {}
    groups_dict = {}
    for f in files:
        k = f.name.split('.')[0]
        #     get index if any controlled vocab in filename
        indices = [i for i, s in enumerate(controlled_vocab) if s.lower() in k.lower()]
        prefix = '_'.join(k.split('_')[:11])
        if len(indices)>0:
            short_name = controlled_vocab[indices[0]].lstrip('_')
            group_name = '%s_%s_%s'%(prefix,version_identifier,short_name)
            short_name = k.split('_')[9]+'_'+short_name # Prepending POL to short_name
        else:
            short_name = k.split('_')[12]
            group_name = '%s_%s_%s'%(prefix,version_identifier,short_name)
            short_name = k.split('_')[9]+'_'+k.split('_')[12]  # Prepending POL to short_name
        file_dict[k] = {'file':f,'group':group_name,'short_name':short_name}
        groups_dict[group_name] = {'items':[],'name':group_name,'short_name':short_name}
    df = pd.DataFrame(file_dict).T
    df.index.name = 'filename'
    df.reset_index(inplace=True)#['group'].unique()
    df.drop(columns=['file'],inplace=True)
    for ug in groups_dict.keys():
        for file_key,file_value in file_dict.items():
            if file_value['group'] == ug:
                groups_dict[ug]['items'].append(file_value['file'])

    #STEP 2: MAKE GROUPS
    groups = []
    for group_key,group_values in groups_dict.items():
        g = metob.Group(name=group_key,items=group_values['items'],short_name=group_values['short_name'])
        groups.append(g)        
        for item in g.items:
            print(g.name,g.short_name,item.name)
        print('')
        
    # STEP 3 STORE GROUPS
    metob.store(groups)
    
    # STEP 4. MAKE SHORTNAMES
    # Make short_filename and short_samplename 
    short_filename_delim_ids = [0,2,4,5,7,9,14]
    short_samplename_delim_ids = [9,12,13,14]
    short_names_df = pd.DataFrame(columns=['sample_treatment','short_filename','short_samplename'])
    ctr = 0
    for f in files:
        short_filename = []
        short_samplename = []
        tokens = f.name.split('.')[0].split('_')
        for id in short_filename_delim_ids:
            short_filename.append(str(tokens[id]))
        for id in short_samplename_delim_ids:
            short_samplename.append(str(tokens[id]))
        short_filename = "_".join(short_filename)
        short_samplename = "_".join(short_samplename)
        short_names_df.loc[ctr, 'full_filename'] = f.name.split('.')[0]
        short_names_df.loc[ctr, 'sample_treatment'] = str(tokens[12]) # delim 12
        short_names_df.loc[ctr, 'short_filename'] = short_filename
        short_names_df.loc[ctr, 'short_samplename'] = short_samplename
        short_names_df.loc[ctr, 'last_modified'] = pd.to_datetime(f.last_modified,unit='s')
        ctr +=1
    short_names_df.sort_values(by='last_modified', inplace=True)
    short_names_df.drop(columns=['last_modified'], inplace=True)
    short_names_df.drop_duplicates(subset=['full_filename'], keep='last', inplace=True)
    short_names_df.set_index('full_filename', inplace=True)
    short_names_df.to_csv(os.path.join(os.path.dirname(output_dir), 'short_names.csv'), sep=',', index=True)
        
    # SETUP QC DIR
    output_data_qc = output_dir
    if not os.path.exists(output_data_qc):
        os.makedirs(output_data_qc)
        
    if(polarity=="POS"):
        exl=['NEG']
        QC_template_filename = pos_templates[1]
    else:
        exl=['POS']
        QC_template_filename = neg_templates[1]
     
    # SELECT GROUPS
    groups = dp.select_groups_for_analysis(name = project_name + my_id,most_recent = True,remove_empty = True,include_list = ['QC'], exclude_list = exl)  
    groups = sorted(groups, key=operator.attrgetter('name'))
    
    # FING CORRESPONDING FILES
    file_df = pd.DataFrame(columns=['file','time','group'])
    for g in groups:
        for f in g.items:
            if hasattr(f, 'acquisition_time'):
                file_df = file_df.append({'file':f, 'time':f.acquisition_time,'group':g}, ignore_index=True)
            else:
                file_df = file_df.append({'file':f, 'time':0,'group':g}, ignore_index=True)

    file_df = file_df.sort_values(by=['time'])
    print("the following QC files were found:")
    for file_data in file_df.iterrows():
        print(file_data[1].file.name)
    
    # SELECT QC ATLAS
    atlases = metob.retrieve('Atlas',name=QC_template_filename,username='vrsingan')
    names = []
    print("the following QC atlases were found:")
    for i,a in enumerate(atlases):
        print(i,a.name,pd.to_datetime(a.last_modified,unit='s'),len(a.compound_identifications))
    print("choosing the most recent atlas:")
    myAtlas = atlases[-1]
    atlas_df = ma_data.make_atlas_df(myAtlas)
    atlas_df['label'] = [cid.name for cid in myAtlas.compound_identifications]
    print(myAtlas.name)
    
    # CREATE METATLAS DATASET FROM FILES AND ATLAS
    print("creating metatlas dataset...")
    all_files = []
    for file_data in file_df.iterrows():
        all_files.append((file_data[1].file,file_data[1].group,atlas_df,myAtlas))
    pool = mp.Pool(processes=min(4, len(all_files)))
    t0 = time.time()
    metatlas_dataset = pool.map(ma_data.get_data_for_atlas_df_and_file, all_files)
    pool.close()
    pool.terminate()
    #If you're code crashes here, make sure to terminate any processes left open.
    print(time.time() - t0)
    
    # MAKE RT PLOTS
    rts_df = dp.make_output_dataframe(input_dataset = metatlas_dataset, fieldname='rt_peak', use_labels=True, output_loc = output_data_qc, summarize=True)
    rts_df.to_csv(os.path.join(output_data_qc,"QC_Measured_RTs.csv"))
    rows = int(math.ceil((rts_df.shape[0]+1)/8))
    cols = 8
    fig = plot.figure()

    gs = gridspec.GridSpec(rows, cols, figure=fig, wspace=1, hspace=2)


    rts_df_copy = rts_df.sort_values(by='standard deviation', ascending=False, na_position='last')

    i = 0
    for line, (index, row) in enumerate(rts_df_copy.iterrows()):
        if not np.isnan(row[:-7]).all():
            ax = fig.add_subplot(gs[i])
            ax.tick_params(direction='out', length=1, pad=0.3, width=0.1, labelsize=0.5)
            ax.scatter(range(rts_df.shape[1]-7),row[:-7], s=0.2)
            i += 1
            ticks_loc = np.arange(0,len(rts_df.columns)-7 , 1.0)
            ax.set_xlim(-0.5,len(rts_df.columns)-7+0.5)
            ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            ax.set_yticks(ax.get_yticks().tolist())
            ax.set_ylim(np.nanmin(row[:-7])-0.12,np.nanmax(row[:-7])+0.12)
            [i.set_linewidth(0.1) for i in ax.spines.values()]
            ax.set_title(row.name, fontsize=0.1)
            ax.set_xlabel('Files', fontsize=1)
            ax.set_ylabel('Actual RTs', fontsize=1)

    plot.savefig(os.path.join(output_data_qc, 'Compound_Atlas_RTs.pdf'), bbox_inches="tight")
    for i,a in enumerate(rts_df.columns):
        print(i, a)
    
    return rts_df, atlas_df, myAtlas

class RT_model:
    def __init__(self, name, coef, intercept):
        self.name = name
        self.coef = coef
        self.intercept=intercept
    
def create_rt_adjustment_model(selected_column,rts_df,atlas_df):
    actual_rts, pred_rts, polyfit_rts = [],[],[]

    current_actual_df = rts_df.loc[:,rts_df.columns[selected_column]]
    bad_qc_compounds = np.where(~np.isnan(current_actual_df))
    current_actual_df = current_actual_df.iloc[bad_qc_compounds]
    current_pred_df = atlas_df.iloc[bad_qc_compounds][['rt_peak']]
    actual_rts.append(current_actual_df.values.tolist())
    pred_rts.append(current_pred_df.values.tolist())

    ransac = RANSACRegressor(random_state=42)
    rt_model_linear = ransac.fit(current_pred_df, current_actual_df)
    coef_linear = rt_model_linear.estimator_.coef_[0]
    intercept_linear = rt_model_linear.estimator_.intercept_

    poly_reg = PolynomialFeatures(degree=2)
    X_poly = poly_reg.fit_transform(current_pred_df)
    rt_model_poly = LinearRegression().fit(X_poly, current_actual_df)
    coef_poly = rt_model_poly.coef_
    intercept_poly = rt_model_poly.intercept_

    for i in range(rts_df.shape[1]-5):
        current_actual_df = rts_df.loc[:,rts_df.columns[i]]
        bad_qc_compounds = np.where(~np.isnan(current_actual_df))
        current_actual_df = current_actual_df.iloc[bad_qc_compounds]
        current_pred_df = atlas_df.iloc[bad_qc_compounds][['rt_peak']]
        actual_rts.append(current_actual_df.values.tolist())
        pred_rts.append(current_pred_df.values.tolist())
        
    x = list(itertools.chain(*pred_rts))
    y = list(itertools.chain(*actual_rts))

    rows = int(math.ceil((rts_df.shape[1]+1)/5))
    cols = 5
    fig = plot.figure(constrained_layout=False)

    gs = gridspec.GridSpec(rows, cols, figure=fig)
    plot.rc('font', size=6)
    plot.rc('axes', labelsize=6)
    plot.rc('xtick', labelsize=3)
    plot.rc('ytick', labelsize=3)


    for i in range(rts_df.shape[1]-5):
        x = list(itertools.chain(*pred_rts[i]))
        y = actual_rts[i]
    
        ax = fig.add_subplot(gs[i])
        ax.scatter(x, y, s=2)
        ax.plot(np.linspace(0, max(x),100), coef_linear*np.linspace(0,max(x),100)+intercept_linear, linewidth=0.5,color='red')
        ax.plot(np.linspace(0, max(x),100), (coef_poly[1]*np.linspace(0,max(x),100))+(coef_poly[2]*(np.linspace(0,max(x),100)**2))+intercept_poly, linewidth=0.5,color='green')
        ax.set_title("File: "+str(i))
        ax.set_xlabel('predicted RTs')
        ax.set_ylabel('actual RTs')
    
    fig_legend = "FileIndex       FileName"
    for i in range(rts_df.shape[1]-5):
        fig_legend = fig_legend+"\n"+str(i)+"        "+rts_df.columns[i]

    fig.tight_layout(pad=0.5)
    plot.text(0,-0.03*rts_df.shape[1], fig_legend, transform=plot.gcf().transFigure)
    plot.savefig(os.path.join(output_data_qc, 'Actual_vs_Predicted_RTs.pdf'), bbox_inches="tight")
    
    qc_df = rts_df[[rts_df.columns[selected_column]]]
    qc_df = qc_df.copy()
    print("Linear Parameters :", coef_linear, intercept_linear)
    print("Polynomial Parameters :", coef_poly,intercept_poly)

    qc_df.columns = ['RT Measured']
    atlas_df.index = qc_df.index
    qc_df['RT Reference'] = atlas_df['rt_peak']
    qc_df['RT Linear Pred'] = qc_df['RT Reference'].apply(lambda rt: coef_linear*rt+intercept_linear)
    qc_df['RT Polynomial Pred'] = qc_df['RT Reference'].apply(lambda rt: (coef_poly[1]*rt)+(coef_poly[2]*(rt**2))+intercept_poly) 
    qc_df['RT Diff Linear'] = qc_df['RT Measured'] - qc_df['RT Linear Pred']
    qc_df['RT Diff Polynomial'] = qc_df['RT Measured'] - qc_df['RT Polynomial Pred']
    qc_df.to_csv(os.path.join(output_data_qc, "RT_Predicted_Model_Comparison.csv"))
    
    output_data_qc = os.path.join(output_dir,"data_qc")
    with open(os.path.join(output_data_qc,'rt_model_linear.txt'), 'w') as f:
        f.write('coef = {}\nintercept = {}\nqc_actual_rts = {}\nqc_predicted_rts = {}'.format(coef_linear, 
                                                                intercept_linear, 
                                                                ', '.join([g.name for g in groups]),
                                                                myAtlas.name))
        f.write('\n'+repr(rt_model_linear.set_params()))
        
    with open(os.path.join(output_data_qc,'rt_model_poly.txt'), 'w') as f:
        f.write('coef = {}\nintercept = {}\nqc_actual_rts = {}\nqc_predicted_rts = {}'.format(coef_poly, 
                                                                intercept_poly, 
                                                                ', '.join([g.name for g in groups]),
                                                                myAtlas.name))
        f.write('\n'+repr(rt_model_poly.set_params()))

    
    print(qc_df)
    
    model_linear = RT_model("linear",coef_linear,intercept_linear)
    model_poly = RT_model("poly",coef_poly,intercept_poly)
    return model_linear, model_poly

def apply_rt_adjustment(model,project_name,my_id,output_dir,pos_atlas_indices = [0,4],neg_atlas_indices = [0,4],save_to_db = True):
    pp=project_name.split("_")
    free_text = my_id.replace("%","")+"_"+pp[0]+"_"+pp[2]+"_"+pp[3]+"_"+pp[4]+"_"+pp[5] # this will be appended to the end of the csv filename exported
    
    output_data_qc = os.path.join(output_dir,"data_qc")
    
    # ATLAS TEMPLATES (alternatively, can be read from a file?)
    pos_templates = ['HILICz150_ANT20190824_TPL_EMA_Unlab_POS',
    'HILICz150_ANT20190824_TPL_QCv3_Unlab_POS',
    'HILICz150_ANT20190824_TPL_ISv5_Unlab_POS',
    'HILICz150_ANT20190824_TPL_ISv5_13C15N_POS',
    'HILICz150_ANT20190824_TPL_IS_LabUnlab2_POS']

    neg_templates = ['HILICz150_ANT20190824_TPL_EMA_Unlab_NEG',
    'HILICz150_ANT20190824_TPL_QCv3_Unlab_NEG',
    'HILICz150_ANT20190824_TPL_ISv5_Unlab_NEG',
    'HILICz150_ANT20190824_TPL_ISv5_13C15N_NEG',
    'HILICz150_ANT20190824_TPL_IS_LabUnlab2_NEG']

    for ix in pos_atlas_indices:
        atlases = metob.retrieve('Atlas',name=pos_templates[ix], username='vrsingan')
        prd_atlas_name = pos_templates[ix].replace('TPL', 'PRD')
        if free_text != '':
            prd_atlas_name = prd_atlas_name+"_"+free_text
        prd_atlas_filename = prd_atlas_name+'.csv'
        myAtlas = atlases[-1]
        PRD_atlas_df = ma_data.make_atlas_df(myAtlas)
        PRD_atlas_df['label'] = [cid.name for cid in myAtlas.compound_identifications]
        if model.name == 'linear':
            PRD_atlas_df['rt_peak'] = PRD_atlas_df['rt_peak'].apply(lambda rt: model.coef*rt+model.intercept)
        else:
            PRD_atlas_df['rt_peak'] = PRD_atlas_df['rt_peak'].apply(lambda rt: (model.coef[1]*rt)+(model.coef[2]*(rt**2))+model.intercept)
        PRD_atlas_df['rt_min'] = PRD_atlas_df['rt_peak'].apply(lambda rt: rt-.5)
        PRD_atlas_df['rt_max'] = PRD_atlas_df['rt_peak'].apply(lambda rt: rt+.5)
    
        PRD_atlas_df.to_csv(os.path.join(output_data_qc, prd_atlas_filename), index=False)
    
        if save_to_db:
            dp.make_atlas_from_spreadsheet(PRD_atlas_df,
                          prd_atlas_name,
                          filetype='dataframe',
                          sheetname='',
                          polarity = 'positive',
                          store=True,
                          mz_tolerance = 12)
        print(prd_atlas_name+" Created!")

    for ix in neg_atlas_indices:
        atlases = metob.retrieve('Atlas',name=neg_templates[ix], username='vrsingan')
        prd_atlas_name = neg_templates[ix].replace('TPL', 'PRD')
        if free_text != '':
            prd_atlas_name = prd_atlas_name+"_"+free_text
        prd_atlas_filename = prd_atlas_name+'.csv'
        myAtlas = atlases[-1]
        PRD_atlas_df = ma_data.make_atlas_df(myAtlas)
        PRD_atlas_df['label'] = [cid.name for cid in myAtlas.compound_identifications]
        if model.name == 'linear':
            PRD_atlas_df['rt_peak'] = PRD_atlas_df['rt_peak'].apply(lambda rt: model.coef*rt+model.intercept)
        else:
            PRD_atlas_df['rt_peak'] = PRD_atlas_df['rt_peak'].apply(lambda rt: (model.coef[1]*rt)+(model.coef[2]*(rt**2))+model.intercept)
        PRD_atlas_df['rt_min'] = PRD_atlas_df['rt_peak'].apply(lambda rt: rt-.5)
        PRD_atlas_df['rt_max'] = PRD_atlas_df['rt_peak'].apply(lambda rt: rt+.5)
    
        PRD_atlas_df.to_csv(os.path.join(output_data_qc, prd_atlas_filename), index=False)
    
        if save_to_db:
            dp.make_atlas_from_spreadsheet(PRD_atlas_df,
                          prd_atlas_name,
                          filetype='dataframe',
                          sheetname='',
                          polarity = 'negative',
                          store=True,
                          mz_tolerance = 12)
    
        print(prd_atlas_name+" Created!")


def write_run_options(project_directory,project_name,analysis_round,polarity):
    import pickle
    import datetime
    run_options_file=os.path.join(project_directory,project_name,"run_options.py")
    timestamp=datetime.datetime.now().strftime('%Y%m%d:%H%M%S')
    run_options = {
        'analysis_round': analysis_round,
        'polarity': polarity,
        'timestamp': timestamp
    }
    with open(run_options_file, 'wb') as out_file:
        pickle.dump(run_options, out_file)

def read_run_options(project_directory,project_name):
    import pickle
    run_options_file=os.path.join(project_directory,project_name,"run_options.py")
    if os.path.exists(run_options_file):
        run_options = pickle.load(open(run_options_file, 'rb'))
        analysis_round = run_options['analysis_round']
        polarity = run_options['polarity']
        print("run mode is set to : "+analysis_round+","+polarity)
        return analysis_round, polarity
    else:
        print("Run mode has not been set yet. Run a block from 2a before proceeding.")
        return None,None

def pick_atlas(project_name,my_id,analysis_round,polarity,username):
    pp=project_name.split("_")
    free_text = my_id.replace("%","")+"_"+pp[0]+"_"+pp[2]+"_"+pp[3]+"_"+pp[4]+"_"+pp[5] # this is the same format string i used in the RT adjuster workbook

    if(analysis_round=="ISTD"):
        atlas_type="IS"
    if(analysis_round=="FINAL"):
        atlas_type="EMA"
    
    atlases = metob.retrieve('Atlas',name='%PRD_'+atlas_type+'_%'+polarity+'_'+free_text+'%',username=username)
    names = []
    for i,a in enumerate(atlases):
        print(i,a.name,pd.to_datetime(a.last_modified,unit='s'))#len(a.compound_identifications)
    return(atlases)

def filter_atlas(my_atlas,metatlas_dataset,project_name,my_id,output_dir,polarity,num_data_points_passing = 5,peak_height_passing = 4e5):
    atlas_df = ma_data.make_atlas_df(my_atlas)
    atlas_df['label'] = [cid.name for cid in my_atlas.compound_identifications]
    print(my_atlas.name)
    metob.to_dataframe([my_atlas])
    
    atlas_df_passing = dp.filter_atlas(atlas_df=atlas_df, input_dataset=metatlas_dataset, num_data_points_passing = num_data_points_passing, peak_height_passing = peak_height_passing)
    print("# Compounds in Atlas: "+str(len(atlas_df)))
    print("# Compounds passing filter: "+str(len(atlas_df_passing)))
    
    if polarity=="POS":
        pol_string="positive"
    if polarity=="NEG":
        pol_string="negative"

    atlas_passing = my_atlas.name+'_filteredby-datapnts'+str(num_data_points_passing)+'-pkht'+str(peak_height_passing)
    myAtlas_passing = dp.make_atlas_from_spreadsheet(atlas_df_passing,
                          atlas_passing,
                          filetype='dataframe',
                          sheetname='',
                          polarity = pol_string,
                          store=True,
                          mz_tolerance = 12)
    atlases = dp.get_metatlas_atlas(name=atlas_passing,do_print = True, most_recent=True)

def make_metatlas_dataset(my_atlas,final_round,project_name,my_id,output_dir,polarity):
    short_names_df = pd.read_csv(os.path.join(os.path.dirname(output_dir), 'short_names.csv'), sep=',', index_col='full_filename')
    if(polarity=="POS"):
        exl=['NEG','QC','Blank']
    else:
        exl=['POS','QC','Blank']
    
    if(os.path.basename(output_dir)=="FINAL"):
        exl.append('InjB')
    
    print(exl)
    groups = dp.select_groups_for_analysis(name = project_name + my_id,    # <- edit text search string here
                                       most_recent = True,
                                       remove_empty = True,
                                       include_list = [], 
                                       exclude_list = exl,
                                       do_print=False)
    print("sorted groups")
    groups = sorted(groups, key=operator.attrgetter('name'))
    for i,a in enumerate(groups):
        print(i, a.name)
        
    atlas_df = ma_data.make_atlas_df(my_atlas)
    atlas_df['label'] = [cid.name for cid in my_atlas.compound_identifications]
    print(my_atlas.name)
    metob.to_dataframe([my_atlas])

    all_files = []
    if(final_round):
        extra_time=0.5
        ma_data.make_data_sources_tables(groups, my_atlas, output_dir, polarity=polarity)
    else:
        extra_time=0.75
    for my_group in groups:
        for my_file in my_group.items:
            #extra_time = 0 
            extra_time = extra_time # .75 for first run, .5 for final
            extra_mz = 0.00
            all_files.append((my_file,my_group,atlas_df,my_atlas,extra_time,extra_mz))
    pool = mp.Pool(processes=min(4, len(all_files)))
    t0 = time.time()
    metatlas_dataset = pool.map(ma_data.get_data_for_atlas_df_and_file, all_files)
    pool.close()
    pool.terminate()
    print(time.time() - t0)
    
    
    return metatlas_dataset


def get_msms_hits(metatlas_dataset,my_atlas,output_dir,ref_loc='/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs_v3.tab'):
      
    import warnings; warnings.simplefilter('ignore')
    
    hits_file=os.path.join(output_dir,my_atlas.name+'_hits.pkl')
    if os.path.exists(hits_file):
        hits = pickle.load(open(hits_file, "rb"))
    else:
        t0 = time.time()
        hits=dp.get_msms_hits(metatlas_dataset,extra_time=True,keep_nonmatches=True,frag_mz_tolerance=.01, ref_loc = ref_loc)
        pickle.dump(hits, open(hits_file, "wb"))
        print(time.time() - t0)
        print('%s%s' % (len(hits),' <- total number of MSMS spectra found in your files'))
        
    return(hits)

def find_isomers(compound_ix,my_atlas):
    atlas_df = ma_data.make_atlas_df(my_atlas)
    atlas_df['label'] = [cid.name for cid in my_atlas.compound_identifications]
    
    cname=atlas_df.iloc[compound_ix].label+" "+atlas_df.iloc[compound_ix].adduct
    cmz=atlas_df.iloc[compound_ix].mz
    crt=atlas_df.iloc[compound_ix].rt_peak
    
    print("*** Target Compound: ***")
    print("    "+str(compound_ix)+","+cname+":    "+str(cmz)+"    "+str(crt))

    print("*** Isomer Compounds: ***")
    for ix,row in atlas_df.iterrows():
        if(row.mz==cmz):
            rowname=row.label+" "+row.adduct
            if not(rowname == cname):
                print("   "+str(ix)+","+row.label+" "+row.adduct+":    "+str(row.mz)+"    "+str(row.rt_peak))

                
def get_compound_notes(compound_ix, polarity, my_atlas, notes_file):
    if os.path.exists(notes_file):
        notes_df=pd.read_csv(notes_file)
    else:
        print(notes_file+" does not exist")
        return(None)
    
    atlas_df = ma_data.make_atlas_df(my_atlas)
    atlas_df['label'] = [cid.name for cid in my_atlas.compound_identifications]
    
    clabel=atlas_df.iloc[compound_ix].label
    print("*** Compound Notes: ***")
    notes=notes_df[(notes_df.label == clabel) & (notes_df.polarity == polarity)]
    return(notes)
  

def remove_marked_compounds(metatlas_dataset,my_atlas,polarity,kept_string="kept"):
    atlas_df = ma_data.make_atlas_df(my_atlas)
    atlas_df['label'] = [cid.name for cid in my_atlas.compound_identifications]
    
    (atlas_all, atlas_kept, atlas_removed) = dp.filter_by_remove(atlas_df, metatlas_dataset)
    print("# Compounds Total: "+str(len(atlas_all)))
    print("# Compounds Kept: "+str(len(atlas_kept)))
    print("# Compounds Removed: "+str(len(atlas_removed)))

    atlasfilename=my_atlas.name+kept_string  # <- enter the name of the atlas to be stored
    if polarity=="POS":
        pol_string="positive"
    if polarity=="NEG":
        pol_string="negative"
    names = dp.make_atlas_from_spreadsheet(atlas_kept, 
                                       atlasfilename,  # <- DO NOT EDIT THIS LINE
                                       filetype='dataframe',
                                       sheetname='',
                                       polarity = pol_string,
                                       store=True,
                                       mz_tolerance = 12
                                      )   
    
def export_results(metatlas_dataset,hits,my_atlas,output_dir,polarity,remove_existing_plots=False):
    
    ## check these folders for existing plots
    plot_dirs=[
    os.path.join(output_dir,polarity+'_boxplot_mz_centroid'),
    os.path.join(output_dir,polarity+'_boxplot_peak_height'),
    os.path.join(output_dir,polarity+'_boxplot_rt_peak'),
    os.path.join(output_dir,polarity+'_compound_EIC_chromatograms'),
    os.path.join(output_dir,polarity+'_data_sheets'),
    os.path.join(output_dir,polarity+'_msms_mirror_plots')
    ]
    
    good_to_go=True
    import shutil
    for my_dir in plot_dirs:
        if os.path.isdir(my_dir):
            dir_files = os.listdir(my_dir)
            if len(dir_files) > 0: 
                if remove_existing_plots:
                    print(my_dir+" is not empty - removing existing files")
                    shutil.rmtree(my_dir)
                else:
                    print(my_dir+" is not empty - remove files yourself or set remove_existing_plots to True")
                    good_to_go=False

    if not(good_to_go):
        print("*** not exporting files ***")
        return
    
    atlas_identifications = dp.export_atlas_to_spreadsheet(my_atlas,os.path.join(output_dir,'%s_%s%s.csv' % (polarity,my_atlas.name,"export")))
    kwargs = {'min_intensity': 1e4,   # strict = 1e5, loose = 1e3
          'rt_tolerance': .5,    #>= shift of median RT across all files for given compound to reference
          'mz_tolerance': 20,      # strict = 5, loose = 25; >= ppm of median mz across all files for given compound relative to reference
          'min_msms_score': .6, 'allow_no_msms': True,     # strict = 0.6, loose = 0.3 <= highest compound dot-product score across all files for given compound relative to reference
          'min_num_frag_matches': 1, 'min_relative_frag_intensity': .001}   # strict = 3 and 0.1, loose = 1, 0.01 number of matching mzs when calculating max_msms_score and ratio of second highest to first highest intensity of matching sample mzs
    scores_df = fa.make_scores_df(metatlas_dataset,hits)
    scores_df['passing'] = fa.test_scores_df(scores_df, **kwargs)

    atlas_df = ma_data.make_atlas_df(my_atlas)
    atlas_df['label'] = [cid.name for cid in my_atlas.compound_identifications]
    pass_atlas_df, fail_atlas_df, pass_dataset, fail_dataset = fa.filter_atlas_and_dataset(scores_df, atlas_df, metatlas_dataset, column='passing')

    fa.make_stats_table(input_dataset = metatlas_dataset, msms_hits = hits, output_loc = output_dir,min_peak_height=1e5,use_labels=True,min_msms_score=0.01,min_num_frag_matches=1,include_lcmsruns = [],exclude_lcmsruns = ['QC'], polarity=polarity)
    scores_df.to_csv(os.path.join(output_dir,'stats_tables',polarity+'_compound_scores.csv'))
    group = 'index' # 'page' or 'index' or None
    save = True
    share_y = True

    short_names_df = pd.read_csv(os.path.join(os.path.dirname(output_dir), 'short_names.csv'), sep=',', index_col='full_filename')
    dp.make_chromatograms(input_dataset=metatlas_dataset, group=group, share_y=share_y, save=save, output_loc=output_dir, short_names_df=short_names_df, short_names_header='short_samplename', polarity=polarity)
    dp.make_identification_figure_v2(input_dataset = metatlas_dataset, msms_hits=hits, use_labels=True, include_lcmsruns = [],exclude_lcmsruns = ['InjBL','InjBl','QC','Blank','blank',"ExCtrl"], output_loc=output_dir,  short_names_df=short_names_df, polarity=polarity)
    
    peak_height = dp.make_output_dataframe(input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_height', output_loc=os.path.join(output_dir,polarity+'_data_sheets'), short_names_df=short_names_df, polarity=polarity, use_labels=True)
    peak_area = dp.make_output_dataframe(input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='peak_area', output_loc=os.path.join(output_dir,polarity+'_data_sheets'), short_names_df=short_names_df, polarity=polarity, use_labels=True)
    mz_peak = dp.make_output_dataframe(input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_peak', output_loc=os.path.join(output_dir,polarity+'_data_sheets'), short_names_df=short_names_df, polarity=polarity, use_labels=True)
    rt_peak = dp.make_output_dataframe(input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [],fieldname='rt_peak', output_loc=os.path.join(output_dir,polarity+'_data_sheets'), short_names_df=short_names_df, polarity=polarity, use_labels=True)
    mz_centroid = dp.make_output_dataframe(input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='mz_centroid', output_loc=os.path.join(output_dir,polarity+'_data_sheets'), short_names_df=short_names_df, polarity=polarity, use_labels=True)
    rt_centroid = dp.make_output_dataframe(input_dataset = metatlas_dataset,include_lcmsruns = [],exclude_lcmsruns = [], fieldname='rt_centroid', output_loc=os.path.join(output_dir,polarity+'_data_sheets'), short_names_df=short_names_df, polarity=polarity, use_labels=True)
    
    dp.make_boxplot_plots(rt_peak, output_loc=os.path.join(output_dir, polarity+'_boxplot_rt_peak'), ylabel="RT Peak")
    dp.make_boxplot_plots(peak_height, output_loc=os.path.join(output_dir, polarity+'_boxplot_peak_height'), ylabel="Peak Height")
    dp.make_boxplot_plots(mz_centroid, output_loc=os.path.join(output_dir, polarity+'_boxplot_mz_centroid'), ylabel="MZ Centroid")
    
    print("*** files exported successfully! ***")
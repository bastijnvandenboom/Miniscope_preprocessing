# Miniscope_preprocessing
Preprocessing pipeline used to preprocess Miniscope data (NoRMCorre and CNMF-E) in the Willuhn lab at the Netherlands Institute for Neuroscience. 

One can run preprocessing for a given session (single mode) or number of sessions consecutive (batch mode). For both options, it runs NoRMCorre on every video, saves non-rigid tiff file, and uses that to run CNMF-e. 

Important scripts
cai_pipeline_bastijn: single mode
cai_pipeline_bastijn_batch: batch mode

Optional scripts
cai_pipeline_bastijn_batch_cnmfe_videos: same as cai_pipeline_bastijn, but in addition saves denoised and demixed video. Run this on a small video file (~1000 frames).
only_cnmf_e_bastijn: only run CNMF-E (use this if CNMF-E keeps crashing during certain session and you want to manually get the right parameters)

This pipeline makes use of other packages (included in this folder):

  rharkes's Fast_Tiff_Write (https://github.com/rharkes/Fast_Tiff_Write)
  flatironinstitute's NoRMCorre (https://github.com/flatironinstitute/NoRMCorre)
  zhoupc's CNMF_E (https://github.com/zhoupc/CNMF_E)

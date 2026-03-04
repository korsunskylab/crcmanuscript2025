###############################################################################
# Imports
###############################################################################
import warnings
warnings.filterwarnings(action = "ignore")
import os
import sys
import shutil
import random
import openslide
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import Polygon
from shapely.validation import make_valid
from matplotlib.patches import Patch
from PIL import Image
from valis import registration, slide_io


Image.MAX_IMAGE_PIXELS = None


os.environ['JAVA_HOME'] = "/usr/libexec/java_home"
os.environ['DYLD_LIBRARY_PATH'] = (
    "/opt/homebrew/Cellar/openslide/4.0.0/lib:" 
    + os.environ.get("DYLD_LIBRARY_PATH", ""))

project_dir = "/path/to/project/dir"
output_dir = os.path.join(project_dir, "outputs", "HALO_pathology_annotation_transfer")
os.makedirs(output_dir, exist_ok=True)


###############################################################################
# Valis Coregistration
###############################################################################

# Define output file paths and create directories if they do not exist
file_dir = f"{project_dir}/outputs/DAPI_HE_coregister/sample_name"
slide_src_dir = f"{file_dir}/slides"
results_dst_dir = f"{file_dir}/registration"
registered_slide_dst_dir = f"{results_dst_dir}/registered_slides"

for folder in [file_dir, slide_src_dir, results_dst_dir, registered_slide_dst_dir]:
    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)


# ---- Load H&E image ---- #
he_level = 1 # This script will default to looking at level 1 of the H&E pyramid
he_dir = f"{project_dir}/data/Motic_H&E" # directory w/ H&E svs files


def get_image_level(image, lidx):
    im = image.read_region(
        (0, 0), 
        image.get_best_level_for_downsample(
            image.level_downsamples[lidx]), 
        image.level_dimensions[lidx])
    return im
image = openslide.OpenSlide(f"{he_dir}/sample_name.svs") 
im = get_image_level(image, he_level)
im.save(f"{slide_src_dir}/sample_name_HE_lvl{he_level}.png") # Save image as png, crop images as necessary to help with alignment


# ---- Run VALIS ---- #
reference_slide = f"{slide_src_dir}/sample_name_HE_lvl{he_level}.png"
registrar = registration.Valis(
    slide_src_dir, results_dst_dir, reference_img_f=reference_slide,
    align_to_reference=True,
    )
rigid_registrar, non_rigid_registrar, error_df = registrar.register()
registration.kill_jvm()

# Check results
# The non-rigid overlap results tend to look better than the rigid overlap results so viewing those
nonrigid_overlap = f"{results_dst_dir}/overlaps/slides_non_rigid_overlap.png"


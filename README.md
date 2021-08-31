Tracking Telescope
=======
A tracking and alignment framework for EUDAQ type files  
Developed for the ETH CMS Pixel Telescope used at PSI

## Installation
Requirements:
- Cmake3.3+
- ROOT6+
```shell
mkdir build
cd build
cmake ..
make
```

## Raw files
The software requires ROOT TTrees with the following branches (1D arrays)
- col
- row
- plane
- adc

Each branch has to have the same number of entries and the "plane" branch indicates from which plane each entry is.

## Usage

General usage:

```shell 
.TrackingTelescope <run_file> <action> <tel_id> <only_tel=0>
```
- run_file: name of the file to align/track
- action: usage option (0, 1 or 2)
- tel_id: telescope ID
- only_tel: track only the telscope (and not the DUT)
  
There are three usage options:

0. tracking: fit the tracks with a straight line based on the alignment and adds tracking data
1. alignment: creates the alignment data with rotations and translations of the planes (stored in [data/alignments.txt](data/alignments.txt))
2. error analyser: calculate the fit uncertainties of each plane

**Each seperate alignment uses a "telescope_id".**

### Create new alignment
Either create new entry in [config/telescopes.txt](config/telescopes.txt) manually or use
```shell
.add_new_telescope <run_path>
```
Each entry looks like this:
```shell 
# TEL NROCS MASK CALIB ZPOS YEAR TYPE COMMENT
 49   4       0   12     2  2019 PAD  # Aug (runs 0-86)
```
- TEL: Telescope iD (each new alignment gets a new telescope ID)
- NROCS: number of planes (automatically extracted in ```.add_new_telescope``` script)
- MASK: number of the mask file ("0" == "0.txt.") in [data/outer_pixel_mask](data/outer_pixel_mask), the default masks all out pixels for a CMS pixel chip
- CALIB: number of the calibration dir ("12" == "telescope12") in  [data/calibrations](data/calibrations) (see next section for more info)
- ZPOS: zposition of the DUT as in [config/z_pos.txt](config/z_pos.txt) (automatically set in the script)
- YEAR: year of the beam test
- TYPE: type of the DUT (PIX, PAD or BCM)  
 
Now run:
```shell 
.TrackingTelescope <run_file> 1 <tel_id> <only_tel=0> (<max_events>) (<n_iter>) (<res_thresh>) (<angle_thresh>) (<sil_dut>)
```
- run_file: name of the file to align
- tel_id: telescope ID
- only_tel: track only the telscope (and not the DUT)

All arguments in paretheses are optional and the default values are taken from [config/align.txt](config/align.txt)

### Add new calibration
- relates to the pulse height calibration of the CMS pixel chips 
- always add new calibration for each telescope ID (recommended)
- move the phCalibration fit files to [data/calibrations/telescope\<ID>](data/calibrations)
- name them ROC0.txt, ROC1.txt, ... in the order the planes appear in the data.
- the is a script to refit the phCalibration_CX.dat files:
```shell 
python/calibrator.py <in_file_name>
```

### Tracking
Once there is an alignment file, the tracking variables can be added to the data.  
Run:
```shell 
.TrackingTelescope <run_file> 0 <tel_id> <only_tel=0> 
```
- creates new ROOT file <run_file>_withTracks.root 
- new file has also the tracking branches
- also creates a lot of plots in plots/<run_number>
- look at plots/<run_number>/index.html 

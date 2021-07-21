# QLA-AstroSat-SXT-Module

The QLA module is meant to generate the Quick look products for the SXT data 
This tool is developed for the SXT payload operation center @ TIFR Mumbai..

Quick Look Analysis and Web-Updater : AstroSat-SXT qla_web_prod_ana_v05.py

CHECK OFFICIAL POC WEBSITE - Quick Look Products

This project is aimed to keep updating the "quick look analysis and web updater module" developed solely for the Soft X-ray telescope payload operation centre (SXT-POC) TIFR Mumbai, serving to the India's first space based multi-wavelength observatory mission "AstroSat". This project started as part of the developer's main project in 2016 while working as a postdoc @ TIFR during 2014-2017. The later updates were supported generously by Center for Space Research, North West University Potchefstroom Campus as part of its contributions to the SXT team. The recent update "v05" released on 17 June 2021 is partly supported by the developer's postdoctoral funding from the NRF at South African Astronomical Observatory (SAAO).
RELEASE NOTE – Soft X-ray Telescope – POC Reference Document

New Updates on Web-based SXT Quick-Look Products Analysis Module

Date: 09 July 2021 20:00 SAST

Current Module : “qla_web_prod_ana_v05.py”

Version No. : v05

Recent Updates On: 17June 2021

Developer: Sunil Chandra

SAAO, Cape Town

NWU, Potchefstroom, South Africa

Contact:

 - Email: sunilc4astrosat@gmail.com chandra@saao.ac.za

 - Phone/Mobile: +91 8451913647 (INDIA)/
                 +27 637912471(SOUTH AFRICA)
Verifications Conducted by (@ POC, TIFR, Mumbai) : - Sunil Chandra - Nilima Kamble - Sandeep L. Vishwakarma

Module to be used by : POC, TIFR, Mumbai, India

Question: Can it be shared outside ? Answer: No, not without proper approval from the developer

Vision of the module: The SXT data streaming from the AstroSat (mission) control center to the payload operation center (POC) at Tata Institute of Fundamental Reserach (TIFR), Mumbai, is verified, cross checked for the data quality, and analysed to generate level-2 events and full-frame science data products. Please refer to the descriptions and references to the published material mentioned on www.tifr.res.in/~astrosat_sxt for more details about the instrument and SXT data. The data quality report (DQR) file and the SXT data (both the level2 & level1 directories) are uploaded on the AstroSat data center at Indian Space Science Data Center (ISSDC) which is then accessible to the proposer/users via AstroBrowse tool. The process of the DQR generation is mediated by mainly two tools A) “sxt_dqr_generator.py” [to write the DQR files using the logfile from “sxtpipeline” run], and B) “qla_web_prod_ana_v05.py” [henceforth QLA module] to generate the quick look products [full-frame raw lightcurve, spectra, and image] and the data structure alongwith the required “.html” files [web-compatible and ready for upload on the POC webpage https://www.tifr.res.in/~astrosat_sxt/HTMLOUTDIR/input.html]. Please note that the plots shared on the above link are obtained from a rough quick look analysis and hence are not useful for any scientific investigations. These are only indicative of the likely quality of the data from the AstroSat observations. This document intends to provide the complete information about “qla_web_prod_ana_v05.py” module and its recent updates. The QLA module performs the following steps :

Reads a number of observation Ids through an input list [an ascii file - OBSLIST]
Reads the QLA configuration through an input “YAML” file
Merges the level-2 events files from different orbits for individual observations listed in OBSLIST
Generates the plots for the full frame integrated lightcurve, spectra, image and variation of the telescope parameters
Generates a fixed data directory structure including the relevant HTML files ready to be uploaded on the webpage
The merged events files, abovementioned science products are also saved for debugging/failsafe run of the QLA module
RELEASE NOTE – Soft X-ray Telescope – POC Reference Document

Highlights of release v05:

The removal of the low energy artefacts from quick look image.
The lightcurve plots are updated in presentations.
The lightcurve ascii data file (“.ipac”), kept as part of the Auxillary Data at POC, now contains MJD as TIME column and not the OFFSET like earlier. Therefore, it enables us to easily combine the lightcurves of a source from all the observations to make a single long-term quick look lightcurve. This can serve as an indicator of the source activity over a period of time and hence will be very helpful in the planning of the future observations.
A flexibility is added to switch between two events merging options [i.e., "Python" or "Julia" tools]. This is easily handled by enabling a particular option in the configualtion file. Note that the Attitude files are not merged by Julia module and hence attitude plot on web shall not be made for this option.
Another flexibility is added to run the module over a list of already merged events files, again parsed through an input ASCII file, similar to the OBSLIST. For this option, the QLA module should be run in development mode [“--devmode=yes” in the command line options]. The list of the merged events files (listing the full path even if those lying in current working directory) is also parsed on command line by “--mrgdlist=MRGDLIST”. This shall be a major benefit for the observations where automated merger script fails to generate merged events.
This major rewriting of the module enables the output data structure directory (in the default “HTMLOUTDIR”) in every steps of the run, not after generating the plots for all the listed observations [as was adopted in v04 and earlier versions of the module]. This was a major drawback of earlier versions, as the entire run was being repeated if any observation ID from the OBSLIST was failing in between. The new update will save a lot of time, and hence ease the efficiency of the POC.
Please note that “--devmode=yes” option is given highest priority and hence the module will not care about OBSLIST inputs in this mode. The module will look for input of the list of already merged events files parsed through “--mrgdlist=MRGDLISTFILE”, if not enabled the module will skip without any warning.
USAGES : 1. For using it in normal [OBSLIST] mode ... leave other inputs as default for most usages

python qla_web_prod_ana_v05.py -c webanalysis.conf -i OBSIDLISTFILE --devmode=no
**2. For using it in development mode ... leave other inputs as default for most usages **

python qla_web_prod_ana_v05.py -c webanalysis.conf --devmode=no – mrgdlist=MRGDLISTFILE
3. For getting the details of the command line options

python qla_web_prod_ana_v05.py -h
=================================================================== Running 'webanalysis' Tool

Task: qla_web_prod_ana_v05.py Version: v05; Release Date: 2021-07-17 Originally Developed By : Dr. Sunil Chandra, with supports from SAAO, Cape Town and NWU Potchefstroom (previously at TIFR Mumbai), Originally developed on 15 june, 2016

Description: This tool is a complete analysis package for SXT quick look products . Note that this product is only meant to be used at POC in TIFR Mumbai. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

===================================================================

Usage: qla_web_prod_ana_v05.py [options]

Options:

-h, --help show this help message and exit

-i OBSIDLIST, --obsidlist=OBSIDLIST   Input File with List of obsIds

-g GRPMIN, --grpmin=GRPMIN    The group min parameter

-a AUTONBIN, --autonbin=AUTONBIN.  The nuber of bins for auto binning the curve

-c CONFIGFILE, --configfl=CONFIGFILE    The directory path ++ for database

-b BINFLAG, --binflag=BINFLAG    The binflag for input 

-e ESTRING, --estring=ESTRING     Channel min cutoff

-m PHAPLTMODE, --phapltmode=PHAPLTMODE Mode for pha plotting Mode if 'plt' it uses PGPLOT otherwise python Default is 'plt'

-r GRDFLAG, --grdflag=GRDFLAG Events Grade Selection

-l LOGFILE, --logname=LOGFILE   The name of output logfile for debugging 

-o FINALOUTTABLE, --fouttbl=FINALOUTTABLE    The IPAC table name for output

-p TASK, --task=TASK    The Task Input for this tools

 --regxymode=REGXYMODE   The XYMODE to be used for XSELECT to extract product

 --imgplttype=IMGPLTTYPE    The pltmode for output image, options are 'ds9' or 'plt[plot][python]'...Default = 'ds9' 

 --devmode=DEVMODE   The devmode for debugging and checking stts of the code

 --mrgdlist=MRGDLIST The list of input merged dirs to be added. Add with. the fullpath in not in working di


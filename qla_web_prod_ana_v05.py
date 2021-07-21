#!/bin/python

'''
>python qla_web_prod_ana_v05.py -h
=========================================================================
	  Running 'webanalysis' Tool 
 	  Task: run_web_ana_v05.py 
 	  Version: v05; Release Date: 2021-07-17 
	  Originally Developed By : 

 	 	 Dr. Sunil Chandra, 
 	 with supports from SAAO, Cape Town and NWU Potchefstroom
 	 	 (previously at TIFR Mumbai), 
 	 Originally developed on 15 june, 2016 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Description: This tool is a complete analysis package for SXT quick look products
.Note that this product is only meant to be used at POC in TIFR Mumbai.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
=========================================================================

Usage: run_web_ana_v05.py [options] 

Options:
  -h, --help            show this help message and exit
  -i OBSIDLIST, --obsidlist=OBSIDLIST
                        Input File with List of obsIds
  -g GRPMIN, --grpmin=GRPMIN
                        The group min parameter
  -a AUTONBIN, --autonbin=AUTONBIN
                        The nuber of bins for auto binning the curve
  -c CONFIGFILE, --configfl=CONFIGFILE
                        The directory path for database
  -b BINFLAG, --binflag=BINFLAG
                        The binflag for input
  -e ESTRING, --estring=ESTRING
                        Channel min cutoff
  -m PHAPLTMODE, --phapltmode=PHAPLTMODE
                        Mode for pha plotting Mode   if 'plt' it uses PGPLOT
                        otherwise python   Default is 'plt'
  -r GRDFLAG, --grdflag=GRDFLAG
                        Events Grade Selection
  -l LOGFILE, --logname=LOGFILE
                        The name of output logfile for debugging
  -o FINALOUTTABLE, --fouttbl=FINALOUTTABLE
                        The IPAC table name for output
  -p TASK, --task=TASK  The Task Input for this tools
  --regxymode=REGXYMODE
                        The XYMODE to be used for XSELECT to extract product
  --imgplttype=IMGPLTTYPE
                        The pltmode for output image, options are 'ds9' or
                        'plt[plot][python]'...Default = 'ds9'
  --devmode=DEVMODE     The devmode for debugging and checking stts of the
                        code
  --mrgdlist=MRGDLIST   The list of input merged dirs to be added. Add with
                        the fullpath in not in working dir..
'''

def print_preamble(inpstr = 'scriptname'):
    import os
    print (inpstr)
    version = (os.path.splitext(scriptname)[0]).split('v')[-1]

    updatereleasedate = '17 July, 2021'

    creationdate = '15 june, 2016'

    preambtxt = "=========================================================================\n"
    preambtxt += "\t  Running 'webanalysis' Tool \n \t  Task: {} \n \t  Version: v{}; Release Date: 2021-07-17 \n".format(scriptname, version)
    preambtxt += "\t  Originally Developed By : \n\n \t \t Dr. Sunil Chandra, \n \t with supports from SAAO, Cape Town and NWU Potchefstroom"
    preambtxt += "\n \t \t (previously at TIFR Mumbai), \n \t Originally developed on {} \n".format(creationdate)
    preambtxt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"
    preambtxt += "Description: This tool is a complete analysis package for SXT quick look products\n." 
    preambtxt += "Note that this product is only meant to be used at POC in TIFR Mumbai.\n"
    preambtxt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    preambtxt += "=========================================================================\n"
    print (preambtxt)


def copyfilefromlist(filelists, srcdir="./", destdir = "destdir") :
    import shutil

    if srcdir != destdir :
        for filename in filelists :
            shutil.copyfile(filename, destdir)
        ststus = 0


def modify( flatImage ):
    #Flat Image
    'Normalize the flat image using the mean'
    mean = flatImage.mean()
    normalizedFlat = flatImage / mean
    
    #Normalize Object Image
    'Normalize the object image using the normalized flat'
    normalizedObject = objectImage / normalizedFlat
    return normalizedObject


def lcplotter(lcf, binsize = 60, en = "0p3to7p0", outfile = "lcout.pdf", LCOutTable = 'lc.ipac') :
	import numpy as np
	import os
	from astropy.table import Table, Column
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from astropy.io import fits
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter

	f = fits.open(lcf, ignore_missing_end=True)
	lcdt = f[1].data
	lchdr = f[0].header
	f.close()

	tstart = lchdr['TSTART']
	tstop = lchdr['TSTOP']
	mjdref = lchdr['MJDREFI']
	srcname = lchdr['object']
	exposure = lchdr['exposure']
	mjdstart = tstart/86400. + mjdref

	mjdref4plt = int(mjdstart)

	enflag = (en.replace("to","-")).replace("p", ".")
	srcname = srcname.replace(" ", "")

	if len(lcdt) > 3 :
		lcdt = lcdt[lcdt['FRACEXP'] >= 0.9]
		lcdtflxm = np.nanmean(lcdt['RATE']); lcdtflsd = np.nanstd(lcdt['RATE'])
		lcdt = lcdt[lcdt['RATE'] >= (lcdtflxm - 3.5*lcdtflsd)]
		lcdt = lcdt[lcdt['RATE'] <= (lcdtflxm + 3.5*lcdtflsd)]

		enstring = enflag
		
		#...writing to file----------------- 
		tblout = Table()
		timeinmjd = (lcdt['TIME'] + tstart)/86400.0 + mjdref
		tblout.add_column(Column(timeinmjd, name = "TIME"), index=0)
		tblout.add_column(Column(lcdt['RATE'], name = "RATE"), index=1)
		tblout.add_column(Column(lcdt['ERROR'], name = "ERROR"), index=2)
		tblout.add_column(Column(lcdt['FRACEXP'], name = "FRACEXP"), index=3)
		tblout['TIME'].unit = 'MET' 
		tblout['RATE'].unit = 'c/s' 
		tblout['ERROR'].unit = 'c/s' 
		tblout.write(LCOutTable, format="ascii.ipac")
		#tblout.close()
		#-----------------------------------

		fig, ax = plt.subplots(1, sharex=True, figsize = [12,9])
		
		#mjdref4plt = int(x[0])
		x = (timeinmjd - mjdref4plt) * 24  #Time in hours
		y = lcdt['RATE']
		ey = lcdt['ERROR']
		print (np.array([x, y])) 
		xlab = "Time [hrs] since MJD {}".format(mjdref4plt)
		ylab = "Rate [c/s]"

		ax.errorbar(x, y, xerr= binsize/(2*60*60), yerr= ey, fmt="*k", markersize = 8, label = "binsize : {} s".format(binsize))

		ylim = [y.min() - 0.3*(y.max() - y.min()),  y.max() + 0.3*(y.max() - y.min())]
		xlim = [x.min() - 0.1*(x.max() - x.min()),  x.max() + 0.2*(x.max() - x.min())]

		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		ax.set_xlabel(xlab, fontsize = 23, fontname = 'Times New Roman')
		ax.set_ylabel(ylab, fontsize = 23, fontname = 'Times New Roman')
		ax.set_title("SXT Lightcurve; Object : {}; E: {} keV; Exp. : {} ks".format(srcname, enflag, round(exposure/1e3,1)), fontsize = 23, fontname = 'Times New Roman')
		#major_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/7.))
		#minor_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/35.))
		#ax.tick_params(axis = 'x', which = 'major', labelsize = 20)
		#ax.tick_params(axis = 'x', which = 'minor', labelsize = 15)
		#ax.set_xticks(np.round(major_ticks,1))
		#ax.set_xticks(minor_ticks, minor = True)
		ax.grid(which = 'major', alpha = 0.4, linewidth=1, color = 'k', linestyle = '--')
		ax.grid(which = 'minor', alpha = 0.3, linewidth=1, color = 'k', linestyle = '-.')
		ax.legend(loc=1, fontsize = 12)
		plt.xticks(fontsize=21, rotation=0)
		plt.yticks(fontsize=21, rotation=0)
		plt.minorticks_on()
		plt.savefig(outfile, bbox_inches='tight')
		fig.clf()
		plt.close()

		if os.path.exists(outfile) :
			status = "LCPlotter succeded in making lighcurve plots with name : \n {}".format(outfile)
		else :
			status = "Warning:: Something fissy happend could not write the curve plot"
			
	else :
		status = "Error:: The length of the lc data is significantly lower (=< 3), so skipping the plotting part"
	#plt.show()
	return status


def qdp2PySpecPlotter(qdpfile = "qdpfile.qdp", outfile="test.png", srcname="TESTONLY", dateobs="2018-12-23", timeobs="23.45", exp=1356, figboxsize = (12,9)):

    import numpy as np
    import matplotlib.pyplot as plt
    import os
 
    phadt = np.loadtxt(qdpfile, skiprows=3)

    fig = plt.figure(figsize = figboxsize)
    ax = fig.add_subplot(1,1,1)
    x = phadt[:,0] ; ex = phadt[:,1]; y = phadt[:,2]; ey = phadt[:,3]
    
    xlab = "Energy [keV]"
    ylab = r"Norm. Counts s$^{-1}$ keV$^{-1}$"

    ax.errorbar(x, y, xerr= ex, yerr= ey, fmt="+k", label = "Date-OBS : {}T{} \n Exp. : {} ks".format(dateobs, timeobs, round(exp/1e3,1)), markersize = 2)
    
    if (y.min()/y.max()) > 1e-3 :
        ylim0 = y.min() - 0.005*(y.max() - y.min())
    else :
        ylim0 = 1e-4

    ylim = [ylim0,  y.max() + 0.4*(y.max() - y.min())]
    xlim = [0.25,  10.0]
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlab, fontsize = 23, fontname = 'Times New Roman', color='black')
    ax.set_ylabel(ylab, fontsize = 23, fontname = 'Times New Roman', color='black')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title("Full Frame SXT Spectrum [an automated product of WebAnalysisTool]", fontsize = 23, fontname = 'Times New Roman', color='black')
    #ax.tick_params(axis = 'x', which = 'major', labelsize = 12)
    #ax.tick_params(axis = 'x', which = 'minor', labelsize = 0)
    #ax.set_xticks([0.20, 0.50, 1.0, 5.0, 10.0])
    #ax.set_xticks(minor_ticks, minor = True)
    #ax.grid(which = 'major', alpha = 0.4, linewidth=1, color = 'k', linestyle = '--')
    #ax.grid(which = 'minor', alpha = 0.3, linewidth=1, color = 'k', linestyle = '-.')
    ax.legend(loc=3, fontsize=12)
    ax.text(0.60, 0.05,'Target Name : {}'.format(srcname), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, bbox=dict(facecolor='none', edgecolor='grey', boxstyle='round, pad=1', linestyle='--'), fontname = 'Times New Roman', fontsize = 15, color = 'black', rotation = 0)

    ax.text(0.350, 0.50,r'$^\dagger$For Quick Look Purpose Only', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, fontname = 'Times New Roman', fontsize = 27, color = 'gainsboro', rotation = 45)

    plt.xticks(fontsize=19, rotation=0)
    plt.yticks(fontsize=19, rotation=0)
    plt.minorticks_on()
    #plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.savefig(outfile, bbox_inches='tight')
    fig.clf()
    plt.close()
    if os.path.exists(outfile) :
        status = "PHAPlotter succeded in making spectral plot with name {}, using qdpPySpecPlotter: \n".format(outfile)
    else :
        status = "Warning :: Something fissy happend, qdpPySpecPlotter could not write the spectral plot\n"

    return status



def run_grppha(inphafile = 'pha.pha', outphafile = 'pha_gr.pha', COMMAND= None, grpmin = 60, backfilename = 'back.pha', 
					respfilename= 'rsp.pha', ancrfilename = 'arf.arf', debug = False) :
						
	import subprocess, os
	import shlex
    
	cmd = 'grppha '
	cmd += " infile={}".format(inphafile) 
	cmd += " outfile={}".format(outphafile)
	cmd += " clobber=true"
	if COMMAND == None :
		COMMAND = "\'group min {} & chkey backfile {} & chkey respfile {} & chkey ancrfile {} & exit\'".format(grpmin, backfilename, respfilename, ancrfilename)
	cmd += " COMM={}".format(COMMAND)
	
	args = shlex.split(cmd)
		
	p = subprocess.Popen(args) 
	output = p.communicate()
	
	if debug: 
		print (cmd)
		print (output)	
		
	if os.path.exists(outphafile) :
		status = 0
	else :
		status = 1
		
	return status  
    

def round_sig(f, p):
    return float(('%.' + str(p) + 'e') % f)
   

def phaplotter(phafname, grpmin = 60, outfile = "phaout.pdf", sxtrspdir = "./", respfilename = "sxt_pc_mat_g0to12.rmf", 
					ancrfilename = "sxt_onaxis_rad05.arf", backfilename = "All6phaadded_spec.pha", 
					pltmode = "plt", pltdev = "PNG", figboxsize = (12,9)) :
						
	import numpy as np
	import os
	from astropy.io import fits
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter

	backfilename = os.path.join(sxtrspdir, backfilename)
	respfilename = os.path.join(sxtrspdir, respfilename)
	ancrfilename = os.path.join(sxtrspdir, ancrfilename)

	grpdphaname = "{}_gr.pha".format(os.path.splitext(phafname)[0])
	
	if os.path.exists(grpdphaname):
		os.remove(grpdphaname)

	try :
		grpphastts = run_grppha(inphafile = phafname, outphafile = grpdphaname, COMMAND= None, grpmin = 60, backfilename = backfilename,
									respfilename= respfilename, ancrfilename = ancrfilename, debug=False)
									
		phaf4plt = grpdphaname
		
	except :
		grpphastts = 999
		phaf4plt = phafname
		
	print ([grpphastts, sxtrspdir])
	phaf = fits.open(phaf4plt, ignore_missing_end=True)
	phafhdr = phaf[0].header
	phafdata = phaf[1].data
	phaf.close()

	tstart = phafhdr['TSTART']
	tstop = phafhdr['TSTOP']
	mjdref = phafhdr['MJDREFI']
	mjdstart = tstart/86400. + mjdref
	dateobs = phafhdr['date-obs']
	timeobs = phafhdr['time-obs']
	srcname = phafhdr['object']
	exp = phafhdr['EXPOSURE']
	# Clearing white space from name & making it uppercase only
	srcname = srcname.replace(" ","")
	timeobs = "{}:{}:{}".format(timeobs.split(":")[0], timeobs.split(":")[1], int(float(timeobs.split(":")[2])))

	exp2prnt = str(round(exp/1e3,1)).split(".")[0]+"p"+str(round(exp/1e3,1)).split(".")[1]+" ks"


	#making tcl script if not in dir...

	if os.path.exists("tcscrpt2run.tcl") :
		os.remove("tcscrpt2run.tcl")
    
	outfileXspec = os.path.splitext(outfile)[0]
    
	#Running xspec for saving the qdp only
	tclf = open("tcscrpt2run.tcl",'w')
	tclf.write("cpd /xw \n")
	tclf.write("setpl e \n")
	tclf.write("data {} \n".format(phaf4plt))
	if grpphastts != 0 :
		tclf.write("back {}\n".format(backfilename))
		tclf.write("response {}\n".format(respfilename))
		tclf.write("arf {}\n".format(ancrfilename))
	tclf.write("ig **-0.3 7.0-** \n")
	tclf.write("ig bad \n")
	tclf.write("plot ldata \n")
	tclf.write("setplot delete all \n")
	tclf.write("setplot device /xw \n")
	tclf.write("setplot command we {} \n".format("tcscrpt2runOut"))
	tclf.write("plot ldata")
	tclf.close()

	# Running tcl scripts for generating qdp for spectrum
	if os.path.exists('tcscrpt2runOut.qdp') :
		os.remove('tcscrpt2runOut.qdp')
	if os.path.exists('tcscrpt2runOut.pco') :
		os.remove('tcscrpt2runOut.pco')

	os.system('xspec tcscrpt2run.tcl')
	phadt = np.loadtxt('tcscrpt2runOut.qdp', skiprows=3)	
	
	if len(np.shape(phadt)) > 1 and len(phadt) >= 3 :
		#print (phadt)
		minNormCnts = np.min(phadt[:,2])
		maxNormCnts = np.max(phadt[:,2])
		pharange = (maxNormCnts - minNormCnts)


		if pltmode == 'plt' :
			if os.path.exists("tcscrpt2run.tcl") :
				os.remove("tcscrpt2run.tcl")
			if os.path.exists('tcscrpt2runOut.qdp') :
				os.remove("tcscrpt2runOut.qdp")
			if os.path.exists('tcscrpt2runOut.pco') :
				os.remove("tcscrpt2runOut.pco")

			tclf = open("tcscrpt2run.tcl",'w')
			tclf.write("cpd /xw \n")
			tclf.write("setpl e \n")
			tclf.write("data {}\n".format(phaf4plt))
			if grpphastts != 0 :
				tclf.write("back {}\n".format(backfilename))
				tclf.write("response {}\n".format(respfilename))
				tclf.write("arf {}\n".format(ancrfilename))
			tclf.write("ig **-0.3 8.0-** \n")
			tclf.write("ig bad \n")
			tclf.write("plot ldata \n")
			tclf.write("setplot device /{} \n".format(pltdev))
			tclf.write("setplot command cpd ./{}/{} \n".format(outfile, pltdev))

			tclf.write('setplot command la t "Full Frame Raw Spectrum : SXT" \n')
			tclf.write('setplot command la f "{}, Date-OBS: {} & Exp.: {} ks" \n'.format(srcname, dateobs.split('T')[0], round(exp/1e3,1) ))
			tclf.write("setplot command cs 1.5 \n")
			tclf.write("setplot command font roman \n")
			tclf.write("setplot command time off \n")
			tclf.write("setplot command rescale x 0.25 11.0  \n")
			tclf.write("setplot command rescale y {} {} \n".format(round_sig(minNormCnts-0.04*pharange,3), round_sig(maxNormCnts+0.09*pharange,4)))
			tclf.write("setplot command vp 0.11 0.12 0.98 0.89 \n")
			tclf.write("setplot command ma 1 on 1 \n")
			tclf.write('setplot command la x "Energy (KeV) " \n')
			tclf.write("plot \n")
			tclf.write("quit  \n")
			tclf.close()

			os.system('xspec tcscrpt2run.tcl')
			if os.path.exists('tcscrpt2runOut.qdp') :
				os.remove('tcscrpt2runOut.qdp')
			if os.path.exists('tcscrpt2runOut.pco') :
				os.remove('tcscrpt2runOut.pco')

			#print ([outfile, os.path.splitext(outfile)[0], os.path.splitext(outfile)[-1]])
            
			if os.path.splitext(outfile)[-1] in [".eps", ".ps"] :
				inplt = outfile
				outfile = "{}.png".format(outfileXspec)
				if os.path.exists(inplt) :
					os.system("convert {} {}".format(inplt, outfile))
					os.remove(inplt)
					#os._exit(0)

			if os.path.exists(outfile) :
				status = "PHAPlotter succeded in making spectral plot with name : \n {}".format(outfile)

			else :
				status = "Warning :: Something fissy happend could not write the spectral plot"

		else :

			if len(phadt) >= 3 :
				outfile2pySpec = "{}.png".format(outfileXspec)
				qdp2PySpecPlotter(qdpfile = "tcscrpt2runOut.qdp", outfile = outfile2pySpec, srcname = srcname, dateobs = dateobs, timeobs = timeobs, exp = exp, figboxsize = figboxsize)
				status = "PhaPlotter Status : 0"
			else :
				status = "Error:: The length of the spectrum files is less than 3, so skipping the plot part..."

	else :
		status = "phaPlotter Status : 999"

	return status, grpdphaname



def FrameProductMaker(evtfile,regdir="./",regfile=None,CHANMIN=30,CHANMAX=701,productstem="SrcNameOut", grade_flag="0-12", curvebinsize = 120, regxymode = "sky"):
    
    import os
    if regfile == None :
        print ("The product is contaminated with corner sources")
        status = "The product is contaminated with corner sources"
        regfile = ''
    else :
        status = "The product is made after excluding the corner sources"

    xco_outstr="{}_scrpt.xco".format(productstem)
    
    regfInp=os.path.join(regdir,regfile)
    
    evtdir="./"
    
    fl_xco=open(str(xco_outstr),'a')
    
    fl_xco.write("xsel\n")
    
    fl_xco.write("read events\n")
    
    fl_xco.write("{}\n".format(evtdir))
    
    fl_xco.write("{}\n".format(evtfile))
    
    fl_xco.write("yes\n")
    
    fl_xco.write("set xyname X Y\n")

    fl_xco.write("set PHANAME PI\n")
    
    fl_xco.write("extract all\n")
    
    fl_xco.write(("filter pha_cutoff {} {} \n").format(CHANMIN,CHANMAX))
    
    fl_xco.write(("filter grade {}\n").format(grade_flag))
    
    #fl_xco.write(("filter region %s \n")%(regfInp))
    
    fl_xco.write("extract image\n")
    
    fl_xco.write("save image {}.img \n".format(productstem))

    if regxymode.upper() in ["RAW", "RAWXY", "RAWX", "RAWY", "RAW_XY", "RAW-XY"] :
        fl_xco.write("set xyname RAWX RAWY\n")
        fl_xco.write("extract all\n")

    fl_xco.write("filter region {} \n".format(regfInp))

    fl_xco.write(("set binsize {} \n").format(curvebinsize))
    
    fl_xco.write("extract curve\n")
    
    fl_xco.write("save curve {}.lc \n".format(productstem))
    
    fl_xco.write("clear pha_cutoff\n")
    
    fl_xco.write("extract spectrum\n")
    
    fl_xco.write("save spectrum {}.pha \n".format(productstem))
    
    fl_xco.write("quit\n")
    
    fl_xco.write("no\n")
    
    fl_xco.close()
    
    os.system('xselect @{}'.format(xco_outstr))
    #os.system("less {}".format(regfInp))
    #os.system("rm -rf {}".format(regfInp))
    return status 


def gkern(kernlen=21, nsig=3):
    import scipy.stats as st
    import numpy as np

    """Returns a 2D Gaussian kernel array."""
    interval = (2*nsig+1.)/(kernlen)
    x = np.linspace(-nsig-interval/2., nsig+interval/2., kernlen+1)
    kern1d = np.diff(st.norm.cdf(x))
    kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
    kernel = kernel_raw/kernel_raw.sum()
    return kernel

def imgaeplotter(infile = "imgname", outfile = "outfile.png", plttype = 'ds9') :
    import matplotlib.pyplot as plt
    import os, seaborn, shlex
    import scipy, subprocess
    import scipy.signal
    import numpy as np
    from astropy.io import fits

    if plttype.upper() in ["PLT", "PYTHON", "PLOT"] :
        status = "The input file : {}, is present ...\n".format(infile)
        imgf = fits.open(infile, ignore_missing_end=True)
        imgfdata = imgf[0].data
        imgfhdr = imgf[0].header
        imgf.close()
    
        mrgd_exp = imgfhdr['EXPOSURE']
        srcname = imgfhdr['object']
    
        if len(srcname.split(" ")) > 1 :
            sname = srcname.split(" ")[0]
            for u in range(len(srcname.split(" "))-1) :
                sname = "{}{}".format(sname,srcname.split(" ")[u+1])
            srcname = sname

        r = scipy.signal.convolve2d(imgfdata, gkern(25), mode='same')
        imgf.close()

        seaborn.set_style("whitegrid")
        fig, ax = plt.subplots(figsize=(10,8))
        sub = ax
        sub.set_title("SXT image: Object : {}; E : 0.3-8.0 keV, Exp: {} ks".format(srcname, str(round(mrgd_exp/1e3,1))), fontsize = 21, fontname = 'Times New Roman')
        pos = sub.imshow(r, interpolation='nearest', vmin=.0001, vmax=r.max()-0.9*r.max() ,cmap=plt.cm.RdBu)#cmap=plt.cm.PuBuGn#.afmhot
        fig.colorbar(pos, ax = sub)
        sub.set_xlim(0,1024)
        sub.set_ylim(0,1024)
        sub.set_ylabel("Y [sxy pixel]", fontsize = 21, fontname = 'Times New Roman')
        sub.set_xlabel("X [sky pixel]", fontsize = 21, fontname = 'Times New Roman')
        seaborn.set_style("ticks")
        print ("The Device has made your plot...save it as per your choice...")
        plt.tick_params(axis='both',which='major',labelsize=12)
        plt.tick_params(axis='both',which='minor',labelsize=8)
        plt.xticks(fontsize=19, rotation=0)
        plt.yticks(fontsize=19, rotation=0)
        plt.minorticks_on()
        plt.savefig("{}".format(outfile), bbox_inches='tight')
        fig.clf()
        plt.close()
    else :
        if os.path.exists(infile) :
            status = "The input file : {}, is present ...\n".format(infile) 

            cmd = "ds9 {} -bin about 3800 3800 -bin factor 2 -scale zscale -cmap AIPS0 -colorbar yes ".format(infile)
            cmd += " -colorbar horizontal -colorbar numerics yes -colorbar fontsize 10 -smooth yes -smooth function gaussian"
            cmd += " -smooth radius 5 -smooth sigma 2.5 -crop 3d 0.02 1.5 wcs -grid yes -grid type publication -grid grid style 1"
            cmd += " -grid axes type exterior -grid grid color white -grid tickmarks color blue -grid numerics color red"
            cmd += " -grid numerics type interior -grid numerics font courier -grid numerics fontsize 9 -grid numerics fontweight bold"
            cmd += " -grid title gap 9 -grid title fontsize 15 -grid title font courier -grid title fontweight bold -grid labels yes"
            cmd += " -grid labels gap1 1 -grid labels gap2 3 -grid labels font courier -grid labels fontsize 17 -grid labels fontweight bold"
            cmd += " -zoom to fit -saveimage {} -exit &".format(outfile)
            
            inpcommand = cmd   

            args = shlex.split(inpcommand)
            #print (args)
            p = subprocess.Popen(args)
            print (p.communicate())
        else :
            print ("The input file : {}, is not present ...".format(infile))
            status = "The input file is not present ...\n"
    if os.path.exists(outfile) :
        status = status + "The image is succesfully written as file: {}\n".format(outfile)
    else :
        status = status + "The image maker is failed \n"
    return status


def time_transform(mjd) :
    from astropy.time import Time
    mjdtime = Time(mjd, format='mjd')
    uttime = mjdtime.isot
    decimalyear = mjdtime.decimalyear
    decimalyear = [round(float(dy),2) for dy in decimalyear if dy]
    return  uttime,  decimalyear

def sxt_met2mjd(timevec) :
    mjdt = []
    for t in timevec :
        mjdt.append((t/86400.) + 55197.0)
    return mjdt

def attplotter(infile = "attfile", outfile = "outfile.png", lctstart = None, lctstop = None) :
    import os, shutil
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits
    attf = fits.open(infile, ignore_missing_end=True)
    attfdata = attf[1].data
    attfhdr  = attf[0].header
    attf.close()
    attm = np.nanmean(attfdata['ang_offset'])
    attstd = np.nanstd(attfdata['ang_offset'])

    attfdata = attfdata[attfdata['ang_offset'] <= (attm + 2.5*attstd)] 
    attfdata = attfdata[attfdata['ang_offset'] >= (attm - 2.5*attstd)]

    if len(attfdata) > 1 :
        attf_time = attfdata['TIME']
        attf_ra = attfdata['roll_ra']
        attf_dec = attfdata['roll_dec']
        attf_offset = 60.*(attfdata['ang_offset'])
        ccd_temp = attfdata['CCD_Temp_1']
        #tstart = attfhdr['tstart']
        #tstop = attfhdr['tstop']
    

        #fig = plt.figure()

        fig, ax = plt.subplots(figsize = [10,8])
        #ax = axarr[0]; ax2 = axarr[1]
        ax.errorbar(attf_time, attf_offset, 0, 0, "*r", markersize=3)
        ax.set_xlabel("Time [MET]", fontsize = 21)
        ax.set_ylabel(r"Ang. Offset; $\theta$ [$arcm$]", fontsize =21)
    
        attf_offset4lim = attf_offset[attf_offset > 0]
        #ylim = [ attf_offset4lim.min() - 0.4*( attf_offset4lim.max() - attf_offset4lim.min() ),  attf_offset4lim.max() + 0.4*(attf_offset4lim.max() - attf_offset4lim.min())]
        ylim = [ np.nanmin(attf_offset) - 0.4*( np.nanmax(attf_offset) - np.nanmin(attf_offset) ),  np.nanmax(attf_offset) + 0.4*(np.nanmax(attf_offset) - np.nanmin(attf_offset))]

        xlim = [np.nanmin(attf_time), np.nanmax(attf_time)]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        major_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/7.))
        minor_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/35.))
        ax.tick_params(axis = 'both', which = 'major', labelsize = 4)
        ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor = True)
        #plt.show()
        ax.axhline(y = np.mean(attf_offset4lim), drawstyle = 'steps', color ='black', lw = 2.5, label = "Mean (M)")
        ax.axhline(y = np.mean(attf_offset4lim) + 2*np.std(attf_offset4lim), linestyle='--',drawstyle = 'steps-post', color ='blue', lw = 2.5, label = r"M $\pm$ 2$\sigma$")
        ax.axhline(y = np.mean(attf_offset4lim) - 2*np.std(attf_offset4lim), linestyle='--', drawstyle = 'steps-post', color ='blue', lw = 2.5)
        ax.axvline(x= lctstart, linestyle = '-.', color = 'green', label = 'Start MET for products', lw = 2.5)
        ax.axvline(x= lctstop, linestyle = '--', color = 'green', label = 'Stop MET for products', lw = 2.5)
    
        ax1 = ax.twiny()
    
        ax1Ticks = ax.get_xticks()
        ax2Ticks = ax1Ticks
        ax1.set_xticks(ax2Ticks)
        ax1.set_xbound(ax.get_xbound())
        ax2ytick = [round(float(x),2) for x in sxt_met2mjd(ax2Ticks) if x]
        ax1.set_xticklabels(ax2ytick)
    
        #print ([ylim[0],ylim[1],int((ylim[1]-ylim[0])/35.)])
        major_ticks = np.arange(ylim[0],ylim[1],(ylim[1]-ylim[0])/7.)
        minor_ticks = np.arange(ylim[0],ylim[1],(ylim[1]-ylim[0])/35.)
        #ax.tick_params(axis = 'Y', which = 'major', labelsize = 4)
        #ax.tick_params(axis = 'Y', which = 'minor', labelsize = 0)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor = True)
        #ax.set_yticks(minor_ticks, minor = True)
        ax.grid(which = 'major', alpha = 0.4, linewidth=1, color = 'k', linestyle = '--')
        ax.grid(which = 'minor', alpha = 0.3, linewidth=1, color = 'k', linestyle = '-.')
        plt.tick_params(labelsize=12)
        ax.legend(loc=2)
        plt.minorticks_on()
        plt.savefig(outfile, bbox_inches='tight')
        #plt.show()
        fig.clf()
        plt.close()
        if os.path.exists(outfile) :
            status = "AttPlotter succeded in making attitude plots with name : \n {}".format(outfile)
        else :
            status = "Warning:: Something fissy happend could not write the attitude plot"
    else :
        status = "Error:: The length of attitude data file not sufficient (=<3) so skipping the plotting part"
    return status



def HndleDtBsnL2dirs(DTBASE = None, obsid = None, logf = None):
	
	import os, glob
	import numpy as np
	
	cwdir = os.getcwd()
	files_obsid = glob.glob(os.path.join(DTBASE,"*"+obsid+"*","*tar.gz"))
	#print (files_obsid)
	logf.write("AstroSat ObsID : {}\n".format(obsid))
           
	# Run only if you have atleast one entry in the file list
	if len(files_obsid) > 0 :
		if str(cwdir) == str(DTBASE) :
			print ("Beaware you are running code in data base directory:\n You may end up with corrupting files")
			logf.write("\t Warning:: ooppsss!!! you are working in the data base directory \n")
			action = raw_input("Enter the action to go or not (Y/N or y/n) : ")
			if action in ['N', 'n', 'not', 'Not', 'NOT'] :
				os._exit()

		else :
			print ("Good that the database directory differs to the current working dirctory")
			logf.write("\t Copying {} tar arcieves to current woring directory\n".format(len(files_obsid)))
			copyfilefromlist(files_obsid, srcdir=DTBASE, destdir = "cwdir")
			logf.write("\t Successfully copied the tar archievs to current dir\n")
                        
	#logf.write("AstroSat ObsID List: {}\n".format(files_obsid))
	for tarfiles in files_obsid :
		os.system("tar xzvf {}".format(tarfiles))
		
	return logf


def TableUpdater(Tbl4Web = None, thisscript = None, HTMLOUTDIR = 'HTMLOUTDIR', htmlresourcedir = 'htmlresourcedir', 
					sunweblink = 'sunweblink', webanalysisscript = 'thisscript', dqr_report = 'dqrflag', dqr_dir = 'dqr_dir', 
					dqr_file_name = 'dqr_file_name', pltmode = 'pltmode', STDBKPFILENAME = 'STDBKPFILENAME', forweblink = 'forweblink', 
					webuplink = 'webuplink') :
	
	from astropy.table import Table, Column, vstack, unique
	import numpy as np
	import os, sys, shutil
	
	cwdir = os.getcwd()
	#Sorting Table in descending order
	Tbl4Web = Tbl4Web[Tbl4Web["MJD_Start"].argsort()][::-1]
	Tbl4Web = unique(Tbl4Web, keys = "ObsID")

	#Subsetting the table for current list of obsid's
	currentTbl = Tbl4Web #subsettingTbl(dttbl=Tbl4Web, filterlist=obsidlist, tblkeyword="ObsID")
	
	#make html for individual sources
	sttsindi = makehtmlfile4indisrc(inptbl = currentTbl, HTMLOUTDIR = HTMLOUTDIR, scdtcpflag = 1, htmlresourcedir = htmlresourcedir, 
										sunweblink = sunweblink, webanalysisscript = thisscript, dqr_report = dqr_report,
										dqr_dir = dqr_dir, dqr_file_name = dqr_file_name, pltmode = pltmode)

	# making the main html table
	sttstable, htmltblstr = html_makerv02(newtbl = Tbl4Web, bkphtmlfile = STDBKPFILENAME, htmlout = "input.html",
											htmlresourcedir = htmlresourcedir, webanalysisscript = thisscript, forweblink = forweblink, 
											webuplink = webuplink, HTMLOUTDIR = HTMLOUTDIR)

	#updating the backup html table
	shutil.copyfile(os.path.join(cwdir, "input.html"), os.path.join(cwdir, STDBKPFILENAME))
	shutil.copyfile(os.path.join(cwdir, "input.html"), os.path.join(HTMLOUTDIR, "input.html"))
	shutil.copyfile(os.path.join(cwdir, STDBKPFILENAME), os.path.join(HTMLOUTDIR, STDBKPFILENAME))
            
	if not os.path.exists(os.path.join(HTMLOUTDIR,"LCoutDir")) :
		os.makedirs(os.path.join(HTMLOUTDIR,"LCoutDir")) 
   
	os.system("mv *lc2plt.ipac {}".format(os.path.join(HTMLOUTDIR,"LCoutDir")))
	

#Main Script...
def main() :
	
	import numpy as np
	import matplotlib.pyplot as plt
	import optparse, os, shutil, glob, sys
	from astropy.table import Table, Column, vstack, unique
	from astropy.io import fits
	import seaborn; seaborn.set()
	import scipy
	import scipy.signal
	import scipy.stats as st
	import yaml

	usage = "usage: %prog [options] "
	parser = optparse.OptionParser(usage)

	parser.add_option("-i", "--obsidlist", dest = "obsidlist", help = "Input File with List of obsIds", default = None)
	parser.add_option("-g", "--grpmin", dest = "grpmin", help = "The group min parameter", default = 15)
	parser.add_option("-a", "--autonbin", dest = "autoNbin", help = "The nuber of bins for auto binning the curve", default = 60)
	parser.add_option("-c", "--configfl", dest = "configfile", help = "The directory path for database", default = 'webanalysis.conf')
	parser.add_option("-b", "--binflag", dest = "binflag", help = "The binflag for input", default = 'fixed')
	parser.add_option("-e", "--estring", dest = "estring", help = "Channel min cutoff", default="0p3to7p0")
	parser.add_option("-m", "--phapltmode", dest = "phapltmode", help = "Mode for pha plotting Mode \n if \'plt\' it uses PGPLOT otherwise python \n Default is \'plt\'", default='plt')
	parser.add_option("-r", "--grdflag", dest = "grdflag", help = "Events Grade Selection", default="0-12")
	parser.add_option("-l", "--logname", dest = "logfile", help = "The name of output logfile for debugging", default = "WebAnaLogFile.log")
	parser.add_option("-o", "--fouttbl", dest = "FinalOutTable", help = "The IPAC table name for output", default = "Tble4Web_FinOut.tbl")
	parser.add_option("-p", "--task", dest = "task", help = "The Task Input for this tools", default = 2)   
	parser.add_option("", "--regxymode", dest = "regxymode", help = "The XYMODE to be used for XSELECT to extract product", default = "raw")
	parser.add_option("", "--imgplttype", dest = "imgplttype", help = "The pltmode for output image, options are 'ds9' or 'plt[plot][python]'...Default = 'ds9' ", default = "ds9")
	parser.add_option("", "--devmode", dest = "devmode", help = "The devmode for debugging and checking stts of the code", default = "yes")
	parser.add_option("", "--mrgdlist", dest = "mrgdlist", help = "The list of input merged dirs to be added. Add with the fullpath in not in working dir..", default = None)
    

	(options, args) = parser.parse_args()
	
	obsidlist, grpmin, autoNbin, configfile = options.obsidlist, options.grpmin, options.autoNbin, options.configfile
	regxymode = options.regxymode
	binflag, estring, grade_flag, taskcode =  options.binflag, options.estring, options.grdflag, int(options.task)
	outlogfile, FinalOutTable, pltmode, imgplttype = options.logfile, options.FinalOutTable, options.phapltmode, options.imgplttype
	devmode, mrgdlist = options.devmode, options.mrgdlist
    
	print (devmode)
	with open(configfile, 'r') as yamlfile :
		configfl = yaml.load(yamlfile, Loader=yaml.FullLoader)

	thisscript = sys.argv[0]
	
	#if str(devmode).upper() in ['YES', 'T', "TRUE"] :
	#	devmode = 'yes'
		
	#Reading crucial parameters from conf file
	DTBASE = configfl["ImpLinks"]["DATABASE"]
	SXTCALDB = configfl["ImpLinks"]["CALDB"]
	MERGERSCRPTSTRING = configfl["ImpLinks"]["MERGERSCRIPT"]
	MERGERSCRIPTTYP = configfl["ImpLinks"]["MERGERSCRIPTTYP"]
	
	SEARCH_STRING = configfl["ImpFiles"]["evtstr"]
	ARFNAME = configfl["ImpFiles"]["arf"]
	RMFNAME = configfl["ImpFiles"]["rmf"]
	BKGNAME = configfl["ImpFiles"]["back"]
	FW_ARFNAME = configfl["ImpFiles"]["fw_arf"] 

	HTMLOUTDIR = configfl["DIRSPECS"]["htmloutdir"]

	STDBKPFILENAME = configfl["STDBKPFILENAME"]["bkphtmlfile"]
	if configfl["DQRINFO"]["dqrflag"] :
		dqrflag = configfl["DQRINFO"]["dqrflag"]
	else : 
		dqrflag = None

	if configfl["DQRINFO"]["dqrdir"] :
		dqr_dir = configfl["DQRINFO"]["dqrdir"]
	else :
		dqr_dir = None

	cssref = configfl["HTMLCONFIG"]["cssref"]
	sunweblink = configfl["HTMLCONFIG"]["sunweblink"]
	htmlresourcedir = configfl["HTMLCONFIG"]["resources"]
	webuplink = configfl["HTMLCONFIG"]["webuplink"]
	forweblink = configfl["HTMLCONFIG"]["weblink_flag"]

	dqr_file_name = None

	if estring != None :
		estring = "0p3to7p0"

	chanmin = int(float(estring.split("to")[0].replace('p','.'))*100) + 1
	chanmax = int(float(estring.split("to")[-1].replace('p','.'))*100) + 1

	#CWD
	cwdir = os.getcwd()
    
	txt2prnt_term = "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
	txt2prnt_term += "\t \t Printing the Default input parameters .....\n"
	txt2prnt_term += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
	txt2prnt_term += "Input File with List of obsIds : {}\n".format(obsidlist)
	txt2prnt_term += "The group min parameter for grppha : {}\n".format(grpmin)
	txt2prnt_term += "The number of bins for auto binning the curve : {}\n".format(autoNbin)
	txt2prnt_term += "The name of input config file; (path of data base)\n".format(configfile)
	txt2prnt_term += "The flag for binning the curve : {}\n".format(binflag)
	txt2prnt_term += "The CHANMIN and CHANMAX for curve : {} & {}\n".format(chanmin, chanmax)
	txt2prnt_term += "The Events grade used for all the products : {}\n".format(grade_flag)
	txt2prnt_term += "The name of output logfile : {}\n".format(outlogfile)
	txt2prnt_term += "The name of the output data table, .ipac format : {}\n".format(FinalOutTable)
	txt2prnt_term += "The full path of the database with level2 tar files : {}\n".format(DTBASE)
	
	print ("Current Working Directory is {} \n".format(cwdir))
	
	if taskcode == 1 :
		txt2prnt_term += "You have choosen to run the html maker script\n Please provided the name of file in IPAC format\n"
        
		tblname = str(raw_input("Please Enter the name of IPAC table with data :\t "))
		tbl = Table.read(tblname)

		# making the main html table
		sttstable, htmltblstr = html_makerv02(newtbl = "NewTable.tbl", bkphtmlfile = "Bkphtml.html", htmlout = "output.html", htmlresourcedir = htmlresourcedir, webanalysisscript = "run_webanalysis_v04.py", forweblink = forweblink, webuplink = webuplink, HTMLOUTDIR = HTMLOUTDIR)

	if taskcode == 2:
		txt2prnt_term += "You have choosen to run the product generator + html maker scripts \n Best Wishes!!!\n"
      
		logf = open(outlogfile, 'a')
    
		logf.write("Date and time of Analysis : \n")  
		
		#print (devmode)
		#Devmode #updated on 15 June 2021
		if str(devmode).upper() in ['YES', 'T', "TRUE"] :
			txt2prnt_term += "The devmode is choosen an option \n"
			#adopting the optional input of the merged events from input list
			
			if mrgdlist != None :
				
				txt2prnt_term += " The use of optional merged input through the list is choosen\n "
				#checking the input merged list file
				
				#Importing mrgdproduct yaml file....
				mrgdprodlistf = yaml.load(open(mrgdlist, "r"), Loader=yaml.FullLoader)
				
				mrgddtdir = mrgdprodlistf['datadir']
				mrgdevtlistf = mrgdprodlistf['data']['MrgdEvts']
				mrgdmkflistf = mrgdprodlistf['data']['MrgdMkfs']
				mrgdmkflistf = mrgdprodlistf['data']['MrgdLbts']
                
				if not type(mrgdmkflistf) in [list, np.ndarray] :
					mrgdmkflistf = [None for k in range(len(mrgdevtlistf))]
					  
				mrgdlbtlistf = mrgdprodlistf['data']['MrgdLbts']
                
				if not type(mrgdlbtlistf) in [list, np.ndarray] :
					mrgdlbtlistf = [None for k in range(len(mrgdevtlistf))]

				
				for mrgdfl, mrgdmkffl, mrgdlbtfl in zip(mrgdevtlistf, mrgdmkflistf, mrgdlbtlistf) :
					
					if mrgdfl :
						
						#Copying events files etc to current working dir....
						if str(mrgddtdir) != str(cwdir) :
							shutil.copyfile(os.path.join(mrgddtdir, mrgdfl), os.path.join(cwdir, mrgdfl))
							
							if mrgdmkffl != None :
								shutil.copyfile(os.path.join(mrgddtdir, mrgdmkffl), os.path.join(cwdir, mrgdmkffl))

							if mrgdlbtfl != None :
								shutil.copyfile(os.path.join(mrgddtdir, mrgdlbtfl), os.path.join(cwdir, mrgdlbtfl))

						fevt = fits.open(mrgdfl)
						datamode = fevt[1].header['datamode']
						objectname = fevt[1].header['object']
						obsid = fevt[1].header['obs_id']
						fevt.close()
						dtmodein = datamode
						txt2prnt_term += " Working on merged events input : {} \n ".format(mrgdfl)
						inpmrgdinpfl = mrgdfl
					    
						pcmrgdprod = [mrgdfl, mrgdmkffl, mrgdlbtfl]
					    
						#using merged evts for events n prod...pc mode
						pctblevtret = run_over_events2(mrgdprod = pcmrgdprod, logf = logf, FinalOutTable = FinalOutTable, binflag = binflag, 
														cwdir = cwdir, estring = estring, grade_flag = grade_flag, grpmin = grpmin,
														configfile = configfile, dtmodein = dtmodein, pltmode = pltmode, regxymode = regxymode,
														imgplttype = imgplttype)
						
						#print(pctblevtret)	
						
						#os._exit(9991)							
						Tbl4Web = pctblevtret
						#update the table
						if not type(Tbl4Web) in [int, float, str] : 								
							TableUpdater(Tbl4Web = Tbl4Web, thisscript = thisscript, HTMLOUTDIR = HTMLOUTDIR, htmlresourcedir = htmlresourcedir, sunweblink = sunweblink, webanalysisscript = thisscript, dqr_report = dqrflag, dqr_dir = dqr_dir,dqr_file_name = dqr_file_name, pltmode = pltmode, STDBKPFILENAME = STDBKPFILENAME, forweblink = forweblink, webuplink = webuplink)
						
						os.system("rm {} *Dt*Exp*.xco *pha".format(mrgdfl))
					
		else :
			
			txt2prnt_term += "The regular mode is choosen an option \n"
            
			#Reading the list file into as a vector...
			obsidlist = np.loadtxt(str(obsidlist),'str')

			loopcnt = 0
			
			for obsid in obsidlist :
            
				logf = HndleDtBsnL2dirs(DTBASE = DTBASE, obsid = obsid, logf = logf)
				
				#PC mode data 
				search_string = SEARCH_STRING
				pcevtlist = glob.glob(os.path.join(cwdir, search_string))
					
				if len(pcevtlist) >= 1 :
					pcevtlist.sort()
					dtmodein = "PC"
					
					#merging the list for pc mode
					pcmrgdprod, pcmrgdprodstts = MergeEvtOut(evtlist = pcevtlist, MERGERSCRPTSTRING = MERGERSCRPTSTRING, 
																MERGERSCRIPTTYP = MERGERSCRIPTTYP, dtmodein = dtmodein, logf = logf)
					print (pcmrgdprod)
					#os._exit(12)
									
					#using merged evts for events n prod...pc mode
					pctblevtret = run_over_events2(mrgdprod = pcmrgdprod, logf = logf, FinalOutTable = FinalOutTable, binflag = binflag, 
														cwdir = cwdir, estring = estring, grade_flag = grade_flag, grpmin = grpmin,
														configfile = configfile, dtmodein = dtmodein, pltmode = pltmode, regxymode = regxymode,
														imgplttype = imgplttype)
					if pctblevtret == 999 :
						print ("The procedure for PC mode over this events file was not succesful")
						pcevtsuccessflag = 0
						
					else :
						pcevtsuccessflag = 1
						
				else :
					pcevtsuccessflag = 0

				
				# For FW mode
				fwsearch_string = search_string.split("PC")[0]+"FW"+search_string.split("PC")[1]
				fwevtlist = glob.glob(os.path.join(cwdir, fwsearch_string))
            
				if len(fwevtlist) >= 1 :
					
					fwevtlist.sort()
					dtmodein = "FW"
					
					#merging the list for fw mode
					fwmrgdprod, fwmrgdprodstts = MergeEvtOut(evtlist = fwevtlist, MERGERSCRPTSTRING = MERGERSCRPTSTRING, 
																MERGERSCRIPTTYP = MERGERSCRIPTTYP, dtmodein = dtmodein)
					#using merged evts for events n prod...fw mode
					fwtblevtret = run_over_events2(mrgdprod = fwmrgdprod, logf = logf, FinalOutTable = FinalOutTable, binflag = binflag, 
														cwdir = cwdir, estring = estring, grade_flag = grade_flag, grpmin = grpmin,
														configfile = configfile, dtmodein = dtmodein, pltmode = pltmode, regxymode = regxymode,
														imgplttype = imgplttype)

					if fwtblevtret == 999 :
						print ("The procedure for FW mode over this events file was not succesful")
						wtevtsuccessflag = 0
						
					else :
						wtevtsuccessflag = 1
						
				else :
					wtevtsuccessflag = 0

				#print ([search_string, cwdir, pcevtlist])
				if pcevtsuccessflag == 1 and wtevtsuccessflag == 1 :
					Tbl4Web = vstack([pctblevtret, fwtblevtret])
					Tbl4Web.sort('MJD_Start')
				elif pcevtsuccessflag == 1 and wtevtsuccessflag == 0 :
					Tbl4Web = pctblevtret
					Tbl4Web.sort('MJD_Start')
				elif pcevtsuccessflag == 0 and wtevtsuccessflag == 1 :
					Tbl4Web = fwtblevtret
					Tbl4Web.sort('MJD_Start')
				else :
					Tbl4Web = Table()
					continue

				#Update the table...		
				TableUpdater(Tbl4Web = Tbl4Web, thisscript = thisscript, HTMLOUTDIR = HTMLOUTDIR, htmlresourcedir = htmlresourcedir, sunweblink = sunweblink, webanalysisscript = thisscript, dqr_report = dqrflag, dqr_dir = dqr_dir, dqr_file_name = dqr_file_name, pltmode = pltmode, STDBKPFILENAME = STDBKPFILENAME, forweblink = forweblink, webuplink = webuplink)
					
				os.system("rm *or*mer*evt *Dt*Exp*.xco *pha")
					
				loopcnt = loopcnt + 1
				
			logf.close()	
	else :
		
		txt2prnt_term += "Wrong Input ...Please use only 1 (html only) or 2 (product generator + html maker)\n"
    
	txt2prnt_term += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"

    #printing above...
	print (txt2prnt_term)

	return 0
        

def subsettingTbl(dttbl=None, filterlist=None, tblkeyword="ObsID"):
    import numpy as np
    from astropy.table import Table, Column, vstack

    if (len(dttbl) >= 2) and (len(filterlist) > 0):
        count = 0
        for filter1 in filterlist :
            tmptbl = dttbl[dttbl[tblkeyword] == filter1]    
            if count == 0:
                outtbl = tmptbl
            else :
                outtbl = vstack([outtbl, tmptbl])
            count = count + 1 
        return outtbl
    else :
        print ("Something Went Wrong in Filtering the Table...")
        return 0


#Not in use, half written stuff
def get_contentFrombkphtml(inphtml = "BkpTable.tbl") :
    import shlex
    num_lines = sum(1 for line in open(inphtml))
    
    inpf = file(inphtml, 'r')
    
    for index, line in enumerate(inpf):
        if len(shlex.split(line)) >= 1:
            if shlex.split(line)[0] == "</thead>" :
                theadindex = str(index)
                print ("+++++++++++++++++++++"+theadindex+"+++++++++++++++++++++")
            if shlex.split(line)[0] == "</table>" :
                tableindex = str(index)
                print ("+++++++++++++++++++++"+tableindex+"+++++++++++++++++++++")

    lpcnt = 0
    inpf = file(inphtml, 'r')
    for index, line in enumerate(inpf):
        print ([index, theadindex, tableindex])
        if (index > theadindex) & (index < tableindex) :
            #print ([index, shlex.split(line)])
            if lpcnt == 0 :
                contentOut = line
            else :
                contentOut += line
            lpcnt = lpcnt + 1

    #print (contentOut)
    #return (contentOut)

#Heart definition module 
def make_reg4sxtmodes(mode = "PC", outregfile = "regfile.reg", regxymode = 'sky') :
    
    import os, sys

    mode = mode.upper()
    regf = open(outregfile,'w')
    regf.write('# Region file format: DS9 version 4.1\n')
    
    if regxymode.upper() == "SKY" :
        
        regf.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        regf.write("image\n")
        XX, YY = 500.99976, 500.99935
    
    elif regxymode.upper() in ["RAWXY", "RAW", "RAW_XY", "RAW-XY"] :
        
        regf.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        regf.write("image\n")
        XX, YY = 300.99976, 300.99935
    
    else :
    
        print ("This part is not implimented yet... This part is supposed to decide the source location automatically...using algorithm developed by Patrick Kilian")
        os._exit(0)

    if mode.upper() == "PC" :
        rad = 290.2329
    else :
        rad = 80.0582

    regf.write("circle({}, {}, {})\n".format(XX, YY, rad))
    regf.close()
    #logf.write("\t\t Region File made ...\n")
    return outregfile


def MergeEvtOut(evtlist = None, MERGERSCRPTSTRING = None, MERGERSCRIPTTYP = "PY", dtmodein = "PC", logf = 'logf'): 
	
	import os, glob, shutil
	import numpy as np
	from astropy.io import fits
	
	mrgstt = "+++++++++++++++++++++++++++++++++++++++++++++\n"
	if evtlist != None :
		
		if MERGERSCRPTSTRING != None :
			
			if len(evtlist) == 0:
				mrgstt += "The eventlist is empty so skipping rest procedure\n"
				logf.write("\t Error:: ***The eventlist is empty so skipping rest procedures for\n{}\n".format(obsid))
				#continue

			else:
				print ("Total of {} events files are found in the list so going for merging these".format(len(evtlist)))
				mrgstt += "Total of {} events files are found in the list so going for merging these".format(len(evtlist))
				
				# Access following only for non-zero length of events list
				DtDirName = evtlist[0].split('/')[-1].split('level2')[0]+'level2'
				# saving the list file as text file
				#print (evtlist)
				np.savetxt('evtlist.lis', evtlist, '%s')
				# check whether written it successfully
				if os.path.exists('evtlist.lis') :
					logf.write("\t The events list file is made successfully \n")
					
				# deciding the source name for naming the out
				fevtf = fits.open(evtlist[0], ignore_missing_end=True)
				srcname = fevtf[0].header['object']
				firstorbit = fevtf[0].header['orb_num']
				obsid2use = fevtf[0].header['OBS_ID']
				fevtf.close()
				
				# Getting rid of white spaces in source name
				srcname = srcname.replace(" ","")
				srcname = srcname.upper()
				
				srcname4files = "{}_{}_{}".format(srcname, dtmodein.upper(), obsid2use)
				logf.write("\t ++++++ The source name being investigated : {}\n".format(srcname))

				# Check Whether Merged events already exist
				#mergedEvtFile = glob.glob("{}*or*mer*_cl.evt".format(srcname4files))

				# merge the events files to make one if number > 1
				if len(evtlist) >= 2 :
                    
					#if len(mergedEvtFile) < 1 :
					if MERGERSCRIPTTYP.upper() in ["PYTHON", "PY", "OLD", "O"] :
						os.system("python {} -f evtlist.lis -s {}".format(MERGERSCRPTSTRING, srcname4files))
						os.system("rm *information.log")
						os.system("rm *VerificationLog.txt")
					elif MERGERSCRIPTTYP.upper() in ["JULIA", "JUL", "JU", "J"] :
						outmrgdevtfname = "{}_or00000_11111_merged_cl.evt".format(srcname4files)
						os.system("python {} -l evtlist.lis -o {}".format(MERGERSCRPTSTRING, outmrgdevtfname))
						
					else :
						print ("Could not understand the type of merger scipt...Exiting")
						os._exit(999)

				# Copy and rename the event file if only one
				else :
					
					modname_evts = "{}or{}_{}_merged_cl.evt".format(srcname4files, firstorbit, firstorbit)
					modname_mkf  = "{}or{}_{}_merged.mkf".format(srcname4files, firstorbit, firstorbit)
					modname_lbt  = "{}or{}_{}_merged.lbt".format(srcname4files, firstorbit, firstorbit)

					shutil.copyfile(evtlist[0], os.path.join(cwdir, modname_evts))
					mkffile2cpy = glob.glob("2*/sxt/*/*.mkf")[0]
					lbtfile2cpy = glob.glob("2*/sxt/*/aux/aux*/*.lbt")[0]
					shutil.copyfile(mkffile2cpy, os.path.join(cwdir, modname_mkf))
					shutil.copyfile(lbtfile2cpy, os.path.join(cwdir, modname_lbt))

				# The following will only be followed if merger has made it succesful...
				os.system("rm -rf 20*.0")
				os.system("rm *evtlist.lis")

				#checking whether the merged/single events file made with proper name
				mergedEvtFile = glob.glob("{}*or*mer*_cl.evt".format(srcname4files))
				mkffilelist = glob.glob("{}_or*_merged.mkf".format(srcname4files))
				lbtfilelist = glob.glob("{}_or*_merged.lbt".format(srcname4files))
                
				if type(mergedEvtFile) in [list, np.ndarray] :
					if len(mergedEvtFile) > 0 :
						mergedEvtFile = mergedEvtFile[0]
					else :
						mergedEvtFile = None
				else :
					mergedEvtFile = None
					 
				if type(mkffilelist) in [list, np.ndarray] :
					if len(mkffilelist) > 0 :
						print (mkffilelist)
						mkffile = mkffilelist[0]
					else :
						mkffile = None
				else :
					mkffile = None
				
				if type(lbtfilelist) in [list, np.ndarray] :
					if len(lbtfilelist) > 0 :
						lbtfile = lbtfilelist[0]
					else :
						lbtfile = None
				else :
					lbtfile = None
					
				prod2ret = [mergedEvtFile, mkffile, lbtfile]
				logf.write("\t The length of list for merged files : {}".format(len(mergedEvtFile)))
				#os._exit()
				mrgdstts2ret = 0
			
		else :
			mrgstt += "The input mergescript string is Nonetype, Nothing to do !!\n"
			prod2ret = [None, None, None]
			mrgdstts2ret = 999

	else :
		mrgstt += "The input evtlist is Nonetype, Nothing to do !!\n"
		prod2ret = [None, None, None]
		mrgdstts2ret = 999
	print (mrgstt)
	
	return [prod2ret, mrgstt]


def run_over_events2(mrgdprod = None, logf = "logf", FinalOutTable = "Table.tbl", 
						binflag = 'auto', cwdir='./', estring = None, grade_flag="0-12", 
						grpmin = 60.0, configfile = "Conf.conf", dtmodein = "PC", 
						pltmode = 'python', regxymode = "sky", imgplttype = 'ds9') :

	from astropy.io import fits
	import numpy as np
	import os, shutil, sys, glob
	from astropy.table import Table, Column, vstack, unique
	import yaml

	# Useful input and conf params
	if estring == None :
		estring = "0p3to7p0"

	chanmin = int(float(estring.split("to")[0].replace('p','.'))*100) + 1
	chanmax = int(float(estring.split("to")[-1].replace('p','.'))*100) + 1

	#print ("---------------Input ESTRING = {} /n Input Channel fitering : CHANMIN - {}\n \t\t CHANMAX = {} \n".format(estring, chanmin, chanmax))
	with open(configfile, "r") as yamlfile :
		configfl = yaml.load(yamlfile, Loader=yaml.FullLoader)

	#Reading crucial parameters from conf file
	DTBASE = configfl["ImpLinks"]["DATABASE"]
	SXTCALDB = configfl["ImpLinks"]["CALDB"]
	MERGERSCRPTSTRING = configfl["ImpLinks"]["MERGERSCRIPT"]
	#SEARCH_STRING = configfl["ImpFiles"]["evtstr"]
		
	if dtmodein.upper() == "PC" :
		ARFNAME = configfl["ImpFiles"]["arf"]
	elif dtmodein.upper() == "FW" :
		ARFNAME = configfl["ImpFiles"]["fw_arf"]
	else :
		print ("Wrong Mode of CCD and hence exiting !!!")
		os._exit()
			
	RMFNAME = configfl["ImpFiles"]["rmf"]
	BKGNAME = configfl["ImpFiles"]["back"]

	
	if mrgdprod != None :
		
		mrgdevtlst, mrgdmkflst, mrgdlbtlst = mrgdprod
		
		if mrgdevtlst :
			#print ("---Test Check Steps ---- {} ".format(mrgdprod))	

			mrgdevtfl, mrgdmkffl, mrgdlbtfl = mrgdevtlst, mrgdmkflst, mrgdlbtlst
			if os.path.exists(mrgdevtfl) :
				
				print (":EXPOSURE Events= {}".format(fits.open(mrgdevtfl)[1].header['exposure']))
				logf.write("\t The merged events file {}\n and Total EXPOSURE = {}".format(mrgdevtfl, fits.open(mrgdevtfl)[1].header['exposure']))
                
				# Reading the header info for output table
				mergedEvtFilef = fits.open(mrgdevtfl, ignore_missing_end=True)
				mergedEvthdr = mergedEvtFilef[0].header
				mrgd_exp = mergedEvthdr['EXPOSURE']
				mrgd_exp2prnt = str(round(mrgd_exp/1e3,1)).split('.')[0]+"p"+str(round(mrgd_exp/1e3,1)).split('.')[1]+"ks"
				mrgd_dateobs = mergedEvthdr['DATE-OBS']
				mrgd_timeobs = mergedEvthdr['TIME-OBS']
				mrgd_timeobsArray = mrgd_timeobs.split(':')
				mrgd_timeobs = "{}:{}:{}".format(mrgd_timeobsArray[0], mrgd_timeobsArray[1], int(float(mrgd_timeobsArray[2])))
				dateobs2tbl = "{}T{}".format(mrgd_dateobs,mrgd_timeobs)
				observer = mergedEvthdr['OBSERVER']
				ra_pnt = mergedEvthdr['RA_PNT']
				dec_pnt = mergedEvthdr['DEC_PNT']
				tstart = mergedEvthdr['TSTART']
				tstop = mergedEvthdr['TSTOP']
				mjdref = mergedEvthdr['MJDREFI']
				mjdstart = tstart/86400. + mjdref
				mjdstop = tstop/86400. + mjdref
				mjdobs = mergedEvthdr['MJD-OBS']
				datamode = mergedEvthdr['DATAMODE']
				mjdrange = "{}-{}".format(round(mjdstart,1), round(mjdstop,1))

				obsid = mergedEvthdr['obs_id']
				
				srcname = (mergedEvthdr['OBJECT'].replace(" ","")).upper()
				obsid2use = mergedEvthdr['OBS_ID']
				
				DtDirName = "AS1{}sxt{}00_level2".format(obsid, datamode)
				
				mergedEvtFilef.close()
				
				srcname4files = "{}_{}_{}".format(srcname, dtmodein.upper(), obsid2use)
				
				#data keeping in tmptbl
                        
				tmptbl = Table()
				tmptbl.add_column(Column([DtDirName], name = 'Data_Folder'), index=0)
				tmptbl.add_column(Column([obsid], name = 'ObsID'), index=1)
				tmptbl.add_column(Column([observer], name = 'Observer'), index=2)
				tmptbl.add_column(Column([srcname], name = 'Object'), index=3)
				tmptbl.add_column(Column([ra_pnt], name = 'RA'), index=4)
				tmptbl.add_column(Column([dec_pnt], name = 'DEC'), index=5)
				tmptbl.add_column(Column([mrgd_exp], name = 'SXT_exposure'), index=6)
				tmptbl.add_column(Column([datamode], name = 'DATAMODE'), index = 7)
				tmptbl.add_column(Column([dateobs2tbl], name = 'Date_Start'), index=8)
				tmptbl.add_column(Column([round(mjdstart,1)], name = 'MJD_Start'), index=9)
				tmptbl.add_column(Column([round(mjdstop,1)], name = 'MJD_Stop'), index=10)
                        
                   
				#if loopcnt == 0 :
				Tbl4Web = tmptbl
				# Reading Old IPAC Table and Loading it ....

				if os.path.exists(os.path.join(cwdir, FinalOutTable)) :
					readoldf = Table.read(os.path.join(cwdir, FinalOutTable), format="ascii.ipac")
					Tbl4Web = vstack((Tbl4Web, readoldf))
					os.remove(os.path.join(cwdir, FinalOutTable))

				Tbl4Web = Tbl4Web[Tbl4Web["MJD_Start"].argsort()][::-1]
				Tbl4Web = unique(Tbl4Web, keys = "ObsID")
				Tbl4Web.write(FinalOutTable,format="ascii.ipac")
                        
				#deciding the binsize for lightcurve output
				timedel = 2.3774729

				if str(binflag) == 'auto' :
					logf.write("\t\t Binsize for curve defined in auto mode\n")
					if mrgd_exp >= 2.*autoNbin :
						binsize = int(mrgd_exp/autoNbin)*timedel
					else :
						binsize = 60*timedel               
				elif str(binflag) == 'fixed' :
					logf.write("\t\t Binsize for curve defined in fixed mode\n")
					if mrgd_exp >= 2*60*timedel :
						binsize = 60*timedel
					else :
						binsize = 30*timedel
				elif str(binflag) == 'interactive' :
					logf.write("\t\t Binsize for curve defined in interactive mode\n Always enter the multiple of 2.377s\n")
					binsize = str(raw_input("Enter the binsize for LC : "))

				else :
					logf.write("\t\t Entry for binsize not given correctly, exiting...\n")
					print ("Enter bin flag Correctly")
					os._exit()

				# Defining the name of output products
				productstem = "{}_Dt{}_Exp{}_AO".format(srcname4files,mrgd_dateobs, mrgd_exp2prnt)
				if os.path.exists('regfile.reg') :
					os.remove('regfile.reg')

				regfile = make_reg4sxtmodes(mode = datamode, outregfile = "regfile.reg", regxymode = regxymode)

				anaStts = FrameProductMaker(mrgdevtfl, regdir="./", regfile = regfile, CHANMIN = chanmin, CHANMAX = chanmax, productstem = productstem, grade_flag = grade_flag, curvebinsize = binsize, regxymode = regxymode)
				logf.write("\t {}\n".format(anaStts))

				#Keeping the copy of raw products in an auxillary directory
				auxsrcspecific = "{}".format(srcname4files)

				if not os.path.exists(auxsrcspecific) :
					os.makedirs(auxsrcspecific)

				# Saving image with region
				imgname = "{}.img".format(productstem)
				imgprod_out = "{}_img.png".format(productstem)

				if os.path.exists(os.path.join(cwdir, imgname)) :
					logf.write("\t The input Auto Image Name is {} \n".format(imgname))
					imgpltStts = imgaeplotter(infile = imgname, outfile = imgprod_out, plttype = imgplttype)
					print (imgprod_out)
					
					logf.write("\t {} \n".format(imgpltStts))
					#logf.close()
					#os._exit(1)
					shutil.move(imgname, os.path.join(auxsrcspecific, imgname)) #moving img to aux dir
                            
				# Saving full frame lightcurve plot
				lcf = "{}.lc".format(productstem)
				if os.path.exists(os.path.join(cwdir, lcf)) :
					lcpltStts = lcplotter(lcf, binsize = binsize, en = estring, outfile = "{}_lc.png".format(productstem), LCOutTable = "{}_lc2plt.ipac".format(productstem))
					logf.write("\t {} \n".format(lcpltStts))
					shutil.move(lcf, os.path.join(auxsrcspecific, lcf)) #moving lc to aux dir

				# Saving full frame spectrum plot
				pltdev = "PNG"
				phafname = "{}.pha".format(productstem)
				specprod_out = "{}_spec.png".format(productstem)			
					
				if os.path.exists(os.path.join(cwdir, phafname)) :
					phapltStts, grpdphaname = phaplotter(phafname, grpmin = grpmin, outfile = specprod_out, sxtrspdir = SXTCALDB,
															respfilename = RMFNAME, ancrfilename = ARFNAME, pltmode = pltmode, 
															pltdev = pltdev, backfilename = BKGNAME)
															
					logf.write("\t {} \n".format(phapltStts))
					shutil.move(phafname, os.path.join(auxsrcspecific, grpdphaname)) #moving pha file to aux dir
					#os._exit(1)
					print (specprod_out)
				shutil.copy(mrgdevtfl, os.path.join(auxsrcspecific, mrgdevtfl)) #events file
                        
				#os.system("rm -rf mergedEvtFile[0]")
				os.system("rm *tcscrpt2runOut.qdp")

				toreturn = Tbl4Web
						
				#----- Compressing the Auxillary Dirs....
				if not os.path.exists("AuxProductsDir") :
					os.makedirs("AuxProductsDir")

				# Saving Attitude plots...
				attprod_out = "{}_att.png".format(productstem)
				if mrgdmkffl != None :
					if os.path.exists(mrgdmkffl) :
						
						if os.path.exists(os.path.join(cwdir, mrgdmkffl)) :
							attpltStts = attplotter(infile = mrgdmkffl, outfile = attprod_out, lctstart = tstart, lctstop = tstop)
							logf.write("\t {} \n".format(attpltStts))
							shutil.copy(mrgdmkffl, os.path.join(auxsrcspecific, mrgdmkffl)) #moving mkf file to aux dir
								
				#saving *lbt files
				if mrgdlbtfl != None :
					
					if os.path.exists(mrgdlbtfl) :
						
						if os.path.exists(os.path.join(cwdir, mrgdlbtfl)) :
							shutil.copy(mrgdlbtfl, os.path.join(auxsrcspecific, mrgdlbtfl)) #moving lbt file to aux dir

				auxsrcspecifictar = "{}.tar.gz".format(auxsrcspecific)
				os.system("tar czvf {} {}".format(auxsrcspecifictar, auxsrcspecific))
				if os.path.exists(os.path.join("AuxProductsDir", auxsrcspecifictar)) :
					os.remove(os.path.join("AuxProductsDir", auxsrcspecifictar))

				shutil.move(auxsrcspecifictar, os.path.join("AuxProductsDir", auxsrcspecifictar))
				shutil.rmtree(auxsrcspecific)

						
			else :
				logf.write("\t The expected merged events file is not found in current dir...check!!!\n")
				toreturn = 999
		
		elif len(mergedEvtFile) == 0:

			logf.write("\t Somehow Merged events file was not generated for following list :...\n {} \n".format(evtlist))
			toreturn = 999
			
		else :

			logf.write("\t Somehow multiple events files for same source found so not running the procedure\n")
			toreturn = 999
			
	return toreturn


def getdatetime() :
    import datetime, getpass
    
    date = datetime.datetime.now()
    username = getpass.getuser().upper()
    if len(str(date.month)) < 2 :
        month = "0{}".format(date.month)
    else :
        month = date.month
    
    if len(str(date.day)) < 2 :
        day = "0{}".format(date.day)
    else :
        day = date.day
    
    if len(str(date.hour)) < 2 :
        hour = "0{}".format(date.hour)
    else :
        hour = date.hour
    
    if len(str(date.minute)) < 2 :
        minute = "0{}".format(date.minute)
    else :
        minute = date.minute
    
    if len(str(date.second)) < 2 :
        second = "0{}".format(date.second)
    else :
        second = date.second
    
    date2print = "{}-{}-{}T{}:{}:{}".format(date.year, month, day, hour, minute, int(second))

    return date2print



#The most recent version of html_maker module, currently in use...

def html_makerv02(newtbl = "NewTable.tbl", bkphtmlfile = "Bkphtml.html", htmlout = "output.html", htmlresourcedir = "../", 
					webanalysisscript = "run_webanalysis_v04.py", forweblink = True, 
					webuplink = "https://www.tifr.res.in/~astrosat_sxt", HTMLOUTDIR = "HTMLOUTDIR") :
    
    #Generate the front-end html table and appending it to backup html table...
 
    from astropy.table import Table, unique
    import datetime, getpass
    import os, glob, shutil
    import numpy as np
    
    #inptbl = Table.read(newtbl, format = 'ascii.ipac')
    
    #print ("The Input updated table ... : {}".format(newtbl))
    print ("The Bkp html file for comparision ...: {}".format(bkphtmlfile))
    print ("The Output html file ...: {}".format(htmlout))

    # Some useful parameters to import the html resources 
    jqueryscript = str(os.path.join(htmlresourcedir, "jquery-1.10.1.min.js"))
    highlightjs = str(os.path.join(htmlresourcedir, "highlight.min.js"))
    tablesortscript = str(os.path.join(htmlresourcedir, "jQuery-Plugin-For-Sortable-Searchable-Tables-Tablesort", "tablesort.js"))

    #Start of the html write up....

    htmlheader = '<!DOCTYPE HTML>\n'
    htmlheader += "<html>\n"
    htmlheader += "<head>\n"
    
    #Linking the table structure, imported from CZTI team's DQR page
    htmlheader += '<! Table structure is adopted from same by CZTI team ...>\n'
    htmlheader += '<link rel="stylesheet" type="text/css" href="{}/astrosat_style.css">\n'.format(htmlresourcedir)

    # Description of the page, Main header of the table
    htmlheader += '<meta name="description" content="Prelim. SXT analysis Reports">\n'
    htmlheader += '<title>Details of Preliminary Analysis by SXT POC </title>\n'
    #htmlheader += '<style>\n'
    #htmlheader += 'body {\n'
    #htmlheader += '\t\t background-image: url("{}/launchofastr.jpg")\n'.format(htmlresourcedir)
    #htmlheader += '}\n'
    #htmlheader += '</style>\n'
    htmlheader += '<meta http-equiv="content-type" content="text/html; charset=UTF-8">\n'
    htmlheader += '<meta name="generator" content="{}">\n'.format(webanalysisscript)
    htmlheader += '<!-- Update date: {} -->\n'.format(getdatetime())

    # Importing the class for the html table, shortable, searchable and countable
    htmlheader += '<table class="table-sort table-sort-search table-sort-show-search-count" style="table-layout: fixed; width: 100%; font-size:0.95em;">\n'

    htmlheader += '<thead>\n\n'

    # Importing the side logos
    htmlheader += '<div class="hleft"><img src="{}/ASTROSAT_LOGO.jpg" width=100%></div>\n'.format(htmlresourcedir)
    htmlheader += '<div class="hright"><img src="{}/first_light1.png" width=100%></div>\n'.format(htmlresourcedir) 

    # Adding few more details of the page and the responsible team members + Few important links
    htmlheader += '<div class="hmid">\n'
    htmlheader += '<thead>\n'
    htmlheader += "<center> \n \t <h1> ASTROSAT SXT </h1> \n \t <h2> Prelim. SXT Data Quality Check Report</h2> \n \t <h2> The quick look data are not to be used for any scientific analysis. These are only indicative of the likely quality of Data </h2> \n\t <p style=\"font-size:16px; color:#538b01; font-weight:bold;\"> Last updated on:<span style=\"color:#FF0000\"> {} </br> <span style=\"color:#538b01\"> Created by:</span> <span style=\"color:#FF0000\"> <a href=https://www.researchgate.net/profile/Sunil_Chandra2 target='_blank' >Sunil Chandra</a> </span></br> <span style=\"color:#538b01;font-style:bold\"> Maintained by:<span style=\"color:#FF0000\"> Nilima Kamble</span></p> \n\n".format(getdatetime()) 

    htmlheader += "<p style=\"font-size:16px; color:#538b01; font-weight:bold;\">Important Links : <a href='https://www.tifr.res.in/~astrosat_sxt/poc.php' target='_blank'> a) SXT POC TIFR</a> <br><br\> <a href='https://astrosat-ssc.iucaa.in/' target='_blank'> b) IUCAA Science Support Cell</a> <br><br\> <a href='https://www.isro.gov.in/astrosat-0' target='_blank'> c) ISRO AstroSat </a></p>\n"

    htmlheader += "</center>\n"
    htmlheader += "</div>\n"

    # Text Style for Column Names in the main table
    textstyle = "style=\"font-size:16px; color:#FF0000; font-weight:bold;\""

    # Adding the Column Names in prescribed format...
    htmlheader += '<tr><th class="table-sort" width=30% {}> Data Folder</th>\n \t <th class="table-sort" width=18% {}>OBSID</th> \n \t <th class="table-sort" width=10% {}>Observer</th> \n \t <th class="table-sort" width=15% {}>Object</th> \n \t <th class="table-sort" width=10% {}>RA</th> \n \t <th class="table-sort" width=10% {}>Dec</th> \n \t <th class="table-sort" width=10% {}>Exposure [s]</th> \n \t <th class="table-sort" width=5% {}>Mode</th> \n \t <th class="table-sort" width=13% {}>Date/Time Start</th> \n \t <th class="table-sort" width=10% {}> MJD Start</th> \n \t <th class="table-sort" width=10% {}> MJD Stop</th>\n'.format(textstyle, textstyle, textstyle, textstyle, textstyle, textstyle, textstyle, textstyle, textstyle, textstyle, textstyle)

    htmlheader += "</tr>\n"
    htmlheader += "</thead>\n\n"
    htmlheader += "<tbody>\n"
    # Running the sub-module to make the rows of table using input IPAC table, correct indexpath should be parsed 
    #Define htmldstdir keyword : destination directory for html output and related files

    htmltblcontent = Tbl2htmlcontent(inptbl = newtbl, HTMLOUTDIR = HTMLOUTDIR, forweblink = forweblink, webuplink = webuplink)

    '''
    if os.path.exists(os.path.join("./",bkphtmlfile)) :
        bkpContent = get_contentFrombkphtml(inphtml = bkphtmlfile)
        htmltblcontent = htmltblcontent + bkpContent
    '''

    # Parsing the Scripts Imported for shaping the table ....
    htmlscript = '<script type="text/javascript" src="{}"></script>\n'.format(jqueryscript)
    htmlscript += '<script src="{}"></script>\n'.format(highlightjs)
    htmlscript += '<script type="text/javascript" src="{}"></script>\n'.format(tablesortscript)
    htmlscript += '<script type="text/javascript">\n'
    #htmlscript += '\t\t // For Demo Purposes\n'
    htmlscript += '\t\t $(function () {\n'
    htmlscript += "\t\t \t $('table.table-sort').tablesort();\n"
    htmlscript += "\t\t \t hljs.initHighlightingOnLoad(); // Syntax Hilighting\n"
    htmlscript += '});\n'
    htmlscript += '</script>\n\n'

    #Adding the second script to enable the 'Top' button for clicking
    htmlscript2 = '<button onclick="topFunction()" id="myBtn" title="Go to top">Top</button>\n'

    # Style for Top button
    htmlscript2 += '<style>\n'
    htmlscript2 += '\t body {\n'
    htmlscript2 += '\t \t  font-family: Arial, Helvetica, sans-serif;\n'
    htmlscript2 += '\t \t  font-size: 15px;\n'
    htmlscript2 += '\t }\n'

    htmlscript2 += '#myBtn {\n'
    htmlscript2 += '\t  display: none;\n'
    htmlscript2 += '\t    position: fixed;\n'
    htmlscript2 += '\t    bottom: 20px;\n'
    htmlscript2 += '\t   right: 30px;\n'
    htmlscript2 += '\t    z-index: 99;\n'
    htmlscript2 += '\t    font-size: 17px;\n'
    htmlscript2 += '\t    border: none;\n'
    htmlscript2 += '\t    outline: none;\n'
    htmlscript2 += '\t    background-color: red;\n'
    htmlscript2 += '\t    color: white;\n'
    htmlscript2 += '\t    cursor: pointer;\n'
    htmlscript2 += '\t    padding: 18px;\n'
    htmlscript2 += '\t    border-radius: 4px;\n'
    htmlscript2 += '}\n'

    htmlscript2 += '#myBtn:hover {\n'
    htmlscript2 += '\tbackground-color: #555;\n'
    htmlscript2 += '\t}\n'

    htmlscript2 += '</style>\n\n'

    #Adding the third script to enable the 'scroll' for the "Top" button
    htmlscript3 = '<script>\n'
    htmlscript3 += '// When the user scrolls down 20px from the top of the document, show the button\n'
    htmlscript3 += 'window.onscroll = function() {scrollFunction()};\n'

    htmlscript3 += 'function scrollFunction() {\n'
    htmlscript3 += '\t    if (document.body.scrollTop > 20 || document.documentElement.scrollTop > 20) {\n'
    htmlscript3 += '\t \t  document.getElementById("myBtn").style.display = "block";\n'
    htmlscript3 += '\t   } else {\n'
    htmlscript3 += '\t \t     document.getElementById("myBtn").style.display = "none";\n'
    htmlscript3 += '\t    }\n'
    htmlscript3 += '}\n'

    htmlscript3 += '// When the user clicks on the button, scroll to the top of the document\n'
    htmlscript3 += '\t function topFunction() {\n'
    htmlscript3 += '\tdocument.body.scrollTop = 0;\n'
    htmlscript3 += '\t document.documentElement.scrollTop = 0;\n'
    htmlscript3 += '\t}\n'
    htmlscript3 += '</script>\n\n'

    # Ending the entries for main table, with start of footer
    htmlfooter = "</tbody>\n"
    htmlfooter += "</table>\n"

    #Formulating the text for complete scripts 1+2+3
    htmlscript = htmlscript + htmlscript2 + htmlscript3
    
    #Adding the script part to the footer
    htmlfooter += "{}".format(htmlscript)

    #Ending the body of html page
    htmlfooter += "</body>\n"

    #Ending the html content
    htmlfooter += "</html>"

    # The complete string with details of html table page
    htmlcomplete = htmlheader + htmltblcontent + htmlfooter

    #Writing the html indexfile 
    htmloutf = open(htmlout, 'w')
    htmloutf.write("{}".format(htmlcomplete))

    if os.path.exists(htmlout) :
        status = 1
    else :
        status = 0

    return status, htmlcomplete


def Tbl2htmlcontent(inptbl = "InpTbl.tbl", HTMLOUTDIR = "HTMLOUTDIR", forweblink = True, 
						webuplink = "https://www.tifr.res.in/~astrosat_sxt/") :
    
    # Generate the table content of html file based on generated ipac table....The header is made by previous module...
    from astropy.table import Table
    import datetime, getpass
    import os, glob, shutil
    import numpy as np

    #inptbl = Table.read(inptbl, format = 'ascii.ipac')
    
    DtDirName = inptbl['Data_Folder']
    ObsId = inptbl['ObsID']
    #DtDirName = ObsId
    observer = inptbl['Observer']
    objectname = inptbl['Object']
    ra = inptbl['RA']
    dec = inptbl['DEC']
    exposure = inptbl['SXT_exposure']
    datestart = inptbl['Date_Start']
    mjd_Start = inptbl['MJD_Start']
    mjd_Stop = inptbl['MJD_Stop']
    datamode = inptbl['DATAMODE']

    
    htmltblcontent = "\n"
    
    HTMLOUTDIRNAME = os.path.basename(os.path.normpath(HTMLOUTDIR))
    print ([HTMLOUTDIR, HTMLOUTDIRNAME])
    #os._exit()
    #Looping over content of table and writing the details as table row...
    for i in range(len(inptbl)) :
        htmldstdtdir = os.path.join(HTMLOUTDIR, DtDirName[i])
        
        webuplink2use = os.path.join(webuplink,  HTMLOUTDIRNAME, DtDirName[i])
        #define a keywork pltdirs : the deeper structures keeping the plots files only under plots subdirectory

        if forweblink in [True, "Yes", "Y", "yes", "y", 1]:
            linkpath = webuplink2use
        else :
            linkpath = htmldstdtdir

        indexpath = os.path.join(linkpath, "index.html")

        htmltblcontent += "<tr><td><a href='{}' target='_blank'>{}</a></td>\n".format(indexpath, DtDirName[i])
        htmltblcontent += "\t <td>{}</td> \n".format(ObsId[i])
        htmltblcontent += "\t <td>{}</td> \n".format(observer[i])
        htmltblcontent += "\t <td>{}</td> \n".format(objectname[i])
        htmltblcontent += "\t <td>{}</td> \n".format(round(ra[i],3))
        htmltblcontent += "\t <td>{}</td> \n".format(round(dec[i],3))
        htmltblcontent += "\t <td>{}</td> \n".format(round(exposure[i],2))
        htmltblcontent += "\t <td>{}</td> \n".format(datamode[i])
        htmltblcontent += "\t <td>{}</td> \n".format(datestart[i])
        htmltblcontent += "\t <td>{}</td> \n".format(mjd_Start[i])
        htmltblcontent += "\t <td>{}</td> \n".format(mjd_Stop[i])
        htmltblcontent += "\t </tr>\n"

    #print (htmltblcontent)
    return htmltblcontent


def makehtmlfile4indisrc(inptbl = "InpTbl.tbl", HTMLOUTDIR = None, scdtcpflag = 1, htmlresourcedir = "https://www.tifr.res.in/~astrosat_sxt/HTMLOUTDIR/htmlresourcedir/", sunweblink = "https://www.researchgate.net/profile/Sunil_Chandra2", webanalysisscript = "run_webanalysis_v04.py", dqr_report = False, dqr_dir = None, dqr_file_name = None, logfile = 'logfile.log', webuplink ="https://www.tifr.res.in/~astrosat_sxt/", forweblink = True, pltmode = "plt" ) :
    from astropy.table import Table
    import datetime, getpass
    import os, glob, shutil
    import numpy as np
    from astropy.io import fits

    cwdir = os.getcwd()

    #DtDirNameList = []
    if HTMLOUTDIR == None :
        HTMLOUTDIR = './'

    if not os.path.exists(HTMLOUTDIR) :
        os.makedirs(HTMLOUTDIR)
    
    #inptbl = Table.read(inptbl, format = 'ascii.ipac')
        
    DtDirName = inptbl['Data_Folder']
    ObsId = inptbl['ObsID']
    #DtDirName = ObsId
    observer = inptbl['Observer']
    objectname = inptbl['Object']
    ra = inptbl['RA']
    dec = inptbl['DEC']
    exposure = inptbl['SXT_exposure']
    datestart = inptbl['Date_Start']
    mjd_Start = inptbl['MJD_Start']
    mjd_Stop = inptbl['MJD_Stop']
    datamode = inptbl['DATAMODE']
    htmltblcontent = "\n"
    
    for i in range(len(inptbl)) :

        #Define htmldstdir keyword : destination directory for html output and related files
        htmldstdtdir = os.path.join(HTMLOUTDIR, DtDirName[i])

        #Check whether htmldstdir exists, if not create it
        if not os.path.exists(htmldstdtdir) :
            os.makedirs(htmldstdtdir)

        #define a keywork pltdirs : the deeper structures keeping the plots files only under plots subdirectory
        pltdirs = os.path.join(htmldstdtdir, "plots")

        #Check if pltdirs exists, if not create it
        if not os.path.exists(pltdirs) :
            os.makedirs(pltdirs)

        #Defining keyword objectname2search used for file output stems..includes object name, obsid and datamode
        objectname2search = "{}_{}_{}".format(objectname[i], datamode[i], ObsId[i])

        # Check Whether the standard products are successfully created or not, if yes copy them to HTML related destinations 
        try :
            imgfilename = glob.glob("{}_*Dt{}_Exp*ks_AO_img.png".format(objectname2search, datestart[i].split('T')[0]))[0]
            shutil.move(imgfilename, os.path.join(pltdirs, imgfilename))
            getimgfile = 1
        except :
            imgfilename = "Img file not found"
        #curve
        try:
            lcfilename = glob.glob("{}_*Dt{}_Exp*ks_AO_lc.png".format(objectname2search, datestart[i].split('T')[0]))[0]
            shutil.move(lcfilename, os.path.join(pltdirs, lcfilename))
            getlcfile = 1
        except : 
            lcfilename = "LC File Not found"

        #spectra
        try :
            specfilename = glob.glob("{}_*Dt{}_Exp*ks_AO_spec.png".format(objectname2search, datestart[i].split('T')[0]))[0]
            shutil.move(specfilename, os.path.join(pltdirs, specfilename))
            getspfile = 1
        except :
            specfilename = "Spec File not found"

        #attitude
        try :
            attfilename = glob.glob("{}_*Dt{}_Exp*ks_AO_att.png".format(objectname2search, datestart[i].split('T')[0]))[0]
            shutil.move(attfilename, os.path.join(pltdirs,attfilename))
            getattfile = 1
        except :
            attfilename = "Att File not found"

        # Some sanity check for script...
        #print ("---------------------",objectname[i],"---------------------")
        #print (objectname2search)
        #print ("{}*_*or*merged_cl.evt".format(objectname2search))
        
        # the corresponding events files in the current working directory
        if len(glob.glob("{}*_or*merged_cl.evt".format(objectname2search))) > 0 :
            mevtfile = glob.glob("{}*_or*merged_cl.evt".format(objectname2search))[0]
            
            if len(glob.glob("{}*_or*merged.mkf".format(objectname2search))) :
                mmkffile = glob.glob("{}*_or*merged.mkf".format(objectname2search))[0]
                #mlbtfile = glob.glob("{}*_or*merged.lbt".format(object[i]))[0]
            else :
                mmkffile = None
        else :
            continue
        # Reading the header or the events file for some crucial information...
        evtf = fits.open(mevtfile, ignore_missing_end=True)
        dateend = evtf[0].header['date-end']
        timeend = evtf[0].header['time-end']
        evtf.close()
        
        # String showing the TIME_END together with DATE_END, to be writting in info part...
        timeend2print = "{}:{}:{}".format(timeend.split(":")[0], timeend.split(":")[1], int(float(timeend.split(":")[2])))
        dateend2print = "{}T{}".format(dateend,timeend2print)

        #Little Details about webanalysis script and its version, will be parsed to the html content in the following steps
        webanalysis_version = list( os.path.splitext(webanalysisscript)[0].split("_v")[-1] )
        webanalysis_version_tmp = "version\t\n0"

        for entry in webanalysis_version :
            webanalysis_version_tmp += ".{}".format(entry) 
        webanalysis_version = webanalysis_version_tmp

        HTMLOUTDIRNAME = HTMLOUTDIR.split("/")[-1]

        webuplink2use = os.path.join(webuplink,  HTMLOUTDIRNAME, DtDirName[i])

        if forweblink in ["Yes", "YES", "Y", True, 1, "yes", "y"] :
            linkpath = os.path.join(webuplink2use, "plots")
        else :
            linkpath = pltdirs

        #Writing the actual index file .....

        indifilehtml = "<html>\n"
        indifilehtml += "<head>\n"

        #Writing the settings for the page and bold title with observation ID
        indifilehtml += '<meta http-equiv="content-type" content="text/html; charset=UTF-8">\n'
        indifilehtml += '<meta name="description" content="ASTROSAT SXT Prelim. Ana.: ObsID {}">\n'.format(ObsId[i])
        #indifilehtml += '<!link rel="stylesheet" type="text/css" href="{}/astrosat_style.css">\n'.format(htmlresourcedir)
        indifilehtml += '<link rel="stylesheet" type="text/css" href="{}/astrosat_style.css">\n'.format(htmlresourcedir)
        indifilehtml += '<!-- Creation date: {} -->\n'.format(getdatetime())
        indifilehtml += '<title>ASTROSAT SXT Prelim. Analysis by POC; ObsID: {}</title>\n'.format(ObsId[i])
        indifilehtml += '</head>\n'

        #Writing the top three/four shortcuts as part of body 
        indifilehtml += '<body>\n'
        indifilehtml += '<div id="nav">\n'
        indifilehtml += '<b> Summary : </b>\n'
        indifilehtml += '<a href="#obsinfo">Observation Details</a> |\n'
        indifilehtml += '<a href="#countratenplots"> Lightcurve + Image</a> |\n'
        indifilehtml += '<a href="#specnattitude">Attitude + Spectrum</a> |\n'
        
        #Adding DQR part in the index file, if DQR related links are properly assigned in config file
        if dqr_report :
            #Check whether DQR files for current observations exist
            if os.path.exists(os.path.join(dqr_dir, dqr_file_name)) :

                #Check whether DQR dir inside HTMLOUTDIR directory structure exists or not
                if not os.path.exists(os.path.join(HTMLOUTDIR, dqr_dir)) :
                    #If not create a blanks directory 
                    os.makedirs(os.path.join(HTMLOUTDIR, dqr_dir))
                #If directory exists check the DQR file exists or not 
                if not os.path.exists(os.path.join(pltdirs, dqr_file_name)) :
                    #In not copy it from DQR_DIR
                    shutil.copyfile(os.path.join(dqr_dir, dqr_file_name), os.path.join(pltdirs, dqr_file_name))

                #dqrlinkstring
                dstdqrdirfile = os.path.join(linkpath, dqr_file_name)

                #writing the link to html file
                indifilehtml += '<a href="{}" target="_blank"> <src="{}" width=100%>Data Quality Report</a> |\n'.format(dstdqrdirfile, dstdqrdirfile)

        #Adding the two side logos  
        indifilehtml += '<br />\n'
        indifilehtml += '</div>\n'

        indifilehtml += '<div class="hleft"><img src="{}/ASTROSAT_LOGO.jpg" width=100%></div>\n'.format(htmlresourcedir)
        indifilehtml += '<div class="hright"><img src="{}/first_light1.png" width=100%></div>\n'.format(htmlresourcedir)

        #Adding some infor as a part of headers, e.g. what is main motive of this page ? and When it was created etc?
        indifilehtml += '<div class="hmid">\n'
        indifilehtml += '<center>\n'
            
        indifilehtml += '\t <a name="top"> </a>\n'
        indifilehtml += '\t <h1>ASTROSAT SXT </h1>\n'
        indifilehtml += '\t <h2>Summary of Prelim. Data Lookup by POC</h2>\n'
        indifilehtml += '\t <h2>ObsID: {}</h2>\n'.format(ObsId[i])

        indifilehtml += "<center> <p style=\"font-size:16px; color:#538b01; font-weight:bold;\"> File Creation Date : <span style=\"color:#FF0000\"> {} </span> </br> <span style=\"color:#538b01\"> The Products are Created Using Web-Based Analysis Tool :</span> <span style=\"color:#FF0000\">{}</span> <span style=\"color:#538b01\">({}) developed by </span> <span style=\"color:#FF0000\"> <a href={} target='_blank' >Sunil Chandra</a> </span></br></p> \n".format(getdatetime(), webanalysisscript, webanalysis_version, sunweblink) 

        indifilehtml += "<p style=\"font-size:16px; color:#538b01; font-weight:bold;\">Important Links : <a href='https://www.tifr.res.in/~astrosat_sxt/poc.php' target='_blank'> a) SXT POC TIFR</a> <br><br\> <a href='https://astrosat-ssc.iucaa.in/' target='_blank'> b) IUCAA Science Support Cell</a> <br><br\> <a href='https://www.isro.gov.in/astrosat-0' target='_blank'> c) ISRO AstroSat </a>\n"

        #indifilehtml += '\t <br/>\n'
        indifilehtml += '\t <hr/>\n'
        indifilehtml += '\t <b><a name=\'obsinfo\'><h2>Observations Details :</h2></a></b>\n'
        indifilehtml += '</center>\n'
        indifilehtml += '</div>\n'
        indifilehtml += '<hr>\n' #<div class='toplink'><a href='#top'>^TOP</a></div>      
        indifilehtml += '<center>\n'
        indifilehtml += '<ul>\n'
        #parsing start/end of observations, obsID, datamode, exposure, target name, observer, RA and Dec as a part of Observations Details
        indifilehtml += '<li><b>date-obs</b>:  {}\n'.format(datestart[i])
        indifilehtml += '</li>\n'
        indifilehtml += '<li><b>date-end</b>:  {}\n'.format(dateend2print)
        indifilehtml += '</li>\n'
        indifilehtml += '<li><b>obs_id</b>:  {}\n'.format(ObsId[i])
        indifilehtml += '</li>\n'
        indifilehtml += '<li><b>Data Mode </b>:  {}\n'.format(datamode[i])
        indifilehtml += '</li>\n'
        indifilehtml += '<li><b>exposure</b>:  {} second\n'.format(exposure[i])
        indifilehtml += '</li>\n'
        indifilehtml += "<li><b>sourceid</b>:  <a href=https://simbad.u-strasbg.fr/simbad/sim-id?protocol=html&Ident={}&NbIdent=1& target='_blank'>{}</a>\n".format(objectname[i].upper(), objectname[i])
        indifilehtml += '</li>\n'
        indifilehtml += '<li><b>observer</b>:  {}\n'.format(observer[i].upper())
        indifilehtml += '</li>\n'
        indifilehtml += '<li><b>ra_pnt</b>:  {}\n'.format(ra[i])
        indifilehtml += '</li>\n'
        indifilehtml += '<li><b>dec_pnt</b>:  {}\n'.format(dec[i])
        indifilehtml += '</li>\n'
        indifilehtml += '</ul>\n'
        indifilehtml += '</center>\n'
    
        #Linking the shortcuts to proper point of the page 
        indifilehtml += '<hr>\n'
        indifilehtml += '<a name="countratenplots"><h3>Lightcurve (count rate v/s time) and Image </h3></a> <div class="toplink"><a href="#top"><h2>^TOP</h2></a></div>\n'
        indifilehtml += '<table border=1>\n'
        indifilehtml += '\t <tr> \n'
        
        indifilehtml += '\t \t <th width=50%>Full Frame Lightcurve</th>\n'
        indifilehtml += '\t \t <th width=50%>Sky Image</th>\n'
        indifilehtml += '\t</tr>\n'


        #Lightcurve
        if os.path.exists(os.path.join(pltdirs, lcfilename)) :
            indifilehtml += '\t <td><a href="{}" target="_blank"><img src="{}" width=100%></a></td>\n'.format(os.path.join(linkpath, lcfilename), os.path.join(linkpath, lcfilename))
        else :
            indifilehtml += '\t <td><a href="" target="_blank"><img src="...Not Made..." width=100%></a></td>\n'

        #Images
        if os.path.exists(os.path.join(pltdirs, imgfilename)) :
            indifilehtml += '\t <td><a href="{}" target="_blank"><img src="{}" width=100%></a></td>\n'.format(os.path.join(linkpath, imgfilename), os.path.join(linkpath, imgfilename))
        else :
            indifilehtml += '\t <td><a href="" target="_blank"><img src="...Not Made..." width=100%></a></td>\n'

        indifilehtml += '</tr>\n'

        #Attitude plot
        if os.path.exists(os.path.join(pltdirs, attfilename)) :
            indifilehtml += '<td><a href="{}" target="_blank"><img src="{}" width=100%></a></td>\n'.format(os.path.join(linkpath, attfilename), os.path.join(linkpath, attfilename))
        else :
            indifilehtml += '\t <td><a href="" target="_blank">Sorry Attitude extraction failed<img src="...Not Made..." width=100%></a></td>\n'

        #Spectrum
        if os.path.exists(os.path.join(cwdir, pltdirs, specfilename)) :
            if pltmode.upper() in ["PLT", "AUTO", "XSPEC"] :
                indifilehtml += '<td><a href="{}" target="_blank"><img src="{}" width=100%></a></td>\n'.format(os.path.join(linkpath, specfilename), os.path.join(linkpath,  specfilename))    
            else :
                indifilehtml += '<td><a href="{}" target="_blank"><img src="{}" width=100%></a></td>\n'.format(os.path.join(linkpath, specfilename), os.path.join(linkpath,  specfilename)) 
        else :
            indifilehtml += '\t <td><a href="" target="_blank"><img src="...Not Made..." width=100%></a></td>\n'
        
        indifilehtml += '</tr>\n'

        #End of Table content
        indifilehtml += '</table>\n'

        #Linking shortcur for spectrum and attitude
        indifilehtml += '<div class="toplink"><a name="specnattitude"></a><a href="#top"><h2>^T</h2></a></div>\n'

        #Writing indifilehtml file (index.html file)
        indihtmlout = open('./index.html', 'w')
        indihtmlout.write("{}".format(indifilehtml))
        indihtmlout.close()
        
        shutil.move(os.path.join('./', 'index.html'), os.path.join(htmldstdtdir,'index.html'))
            
    #print (indifilehtml)
    return 0



if __name__ == "__main__":
	
	import sys
	
	sys.settrace

	scriptname = sys.argv[0]

	print_preamble(inpstr = scriptname)

	sys.setrecursionlimit(500000)

	main()

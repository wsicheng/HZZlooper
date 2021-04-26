import os
import numpy as np
import ROOT
import pyRootPlotMaker as ppm
import MT2PlotUtils as utils

def MT2PlotMaker(rootdir, samples, data, dirname, plots, output_dir=".", exts=["pdf"], tag="", scaleMC=True, suffix = '',
                 datatitle = 'Observed', multisuf=None, systset=None, lumi=None, ratioType=0, signame=None, sig_points=None, gencats=None,
                 ratioRange=None, doComposition=False):
    '''
    rootdir contains output of MT2Looper, samples are names of the .root files,
    data is the name of the data file, dirname is the directory within the root file
    to extract plots from, plots are a list of plot definitions from MT2PlotDefs
    exts is a list of file extensions to produce
    note that dirname can be a '+' separated string of directories to add together
    plots is a vector of sets (hname, logy, xragne, yrange, rbin)
    '''

    h_bkg_vecs = [[] for x in plots]
    h_data = []

    dirnames = [s.strip() for s in dirname.split("+")]
    if type(multisuf) == list:
        dirnames = [dirname+s for s in multisuf]
    # dirname = dirname.replace('+','')

    dogencat = False
    ncats = len(samples)
    if type(gencats) == list and ncats == 1:
        ncats = len(gencats)
        dogencat = True

    if lumi==None and data != None:
        if '16'   in data or 'plots16' in output_dir: lumi = 35.9
        if '17'   in data or 'plots17' in output_dir: lumi = 41.5
        if '18'   in data or 'plots18' in output_dir: lumi = 59.7
        if 'run2' in data or 'plotsn2' in output_dir: lumi = 137.2
    else:
        lumi = 35.9
        if '2016' in rootdir: lumi = 35.9
        if '2017' in rootdir: lumi = 41.5
        if '2018' in rootdir: lumi = 59.7
        if 'run2' in rootdir: lumi = 137.2

    ## deal with suffixes
    if suffix != '' and suffix[0] != '_':
        suffix = '_'+suffix

    systs = [None for x in plots]
    drawSystematicBand = False
    if systset is not None and type(systset)==list:
        drawSystematicBand = True
        systs = [[] for x in plots]
        h_bkg_vecs_syst_up = [[[] for s in systset] for p in plots ]
        h_bkg_vecs_syst_dn = [[[] for s in systset] for p in plots ]

    ## get background histograms
    for icat in range(ncats):
        # get the root file for the given sample. This assumes that frag/fake photons come from qcd_ht.root
        if dogencat:
            fn = os.path.join(rootdir, samples[0]+".root")
        else:
            fn = os.path.join(rootdir, samples[icat]+".root")
        fid = ROOT.TFile(fn)

        for iplot in range(len(plots)):
            vn = plots[iplot][0]
            if dogencat: vn += '_'+gencats[icat]
            if suffix != None: vn += suffix
            vni = vn+multisuf[0] if type(multisuf)==list else vn
            h = fid.Get("{0}/h_{1}".format(dirnames[0],vni))
            if not h: print( h, dirnames[0], vni, fn )

            h_bkg_vecs[iplot].append( h )

            # histogram won't exist if there are no events. Replace it with None, handle later
            if type(h_bkg_vecs[iplot][-1])==type(ROOT.TObject()):
                h_bkg_vecs[iplot][-1] = None
            else:
                h_bkg_vecs[iplot][-1].SetDirectory(0)
                # handle the case with more than one directory
                for idir in range(1, len(dirnames)):
                    vni = vn+multisuf[idir] if type(multisuf)==list else vn
                    h = fid.Get("{0}/h_{1}".format(dirnames[idir],vni))
                    h_bkg_vecs[iplot][-1].Add(h)

                if len(plots[iplot]) >= 5:
                    h_bkg_vecs[iplot][-1].Rebin(plots[iplot][4])

            if drawSystematicBand:
                for isys, sys in enumerate(systset):
                    vnup = "{}_{}Up{}".format(plots[iplot][0], sys, suffix) 
                    vndn = "{}_{}Dn{}".format(plots[iplot][0], sys, suffix) 
                    hsysup = fid.Get("{}/h_{}".format(dirnames[0],vnup))
                    hsysdn = fid.Get("{}/h_{}".format(dirnames[0],vndn))
                    if not hsysup or not hsysdn:
                        h_bkg_vecs_syst_up[iplot][isys].append( None )
                        h_bkg_vecs_syst_dn[iplot][isys].append( None )
                        continue
                    h_bkg_vecs_syst_up[iplot][isys].append( hsysup )
                    h_bkg_vecs_syst_dn[iplot][isys].append( hsysdn )
                    h_bkg_vecs_syst_up[iplot][isys][-1].SetDirectory(0)
                    h_bkg_vecs_syst_dn[iplot][isys][-1].SetDirectory(0)

                    if len(plots[iplot]) >= 5:
                        h_bkg_vecs_syst_up[iplot][isys][-1].Rebin(plots[iplot][4])
                        h_bkg_vecs_syst_dn[iplot][isys][-1].Rebin(plots[iplot][4])

        fid.Close()

    # deal with nonexistent histograms <--
    skipList = []
    for i in range(len(plots)):
        firstGood = -1
        for icat in range(ncats):
            if h_bkg_vecs[i][icat] != None:
                firstGood = icat
                break
        if firstGood==-1:
            # raise RuntimeError("all background histograms are empty for {0}/h_{1}!".format(dirname,plots[i][0]))
            print( "Error: all background histograms are empty for {0}/h_{1}! Skipping!".format(dirname,plots[i][0]) )
            skipList.append(i)
            continue
        for icat in range(ncats):
            if h_bkg_vecs[i][icat] == None:
                h_bkg_vecs[i][icat] = h_bkg_vecs[i][firstGood].Clone()
                h_bkg_vecs[i][icat].Reset()
        if drawSystematicBand:
            h_bkg_tot = h_bkg_vecs[i][firstGood].Clone()
            h_bkg_tot.GetNbinsX()
            syst_up = [0.0] * h_bkg_tot.GetNbinsX()
            syst_dn = [0.0] * h_bkg_tot.GetNbinsX()
            for icat in range(ncats):
                if icat != firstGood and h_bkg_vecs[i][icat] != None:
                    h_bkg_tot.Add(h_bkg_vecs[i][icat])
            for isys, sys in enumerate(systset):
                h_bkg_tot_syst_up = h_bkg_vecs_syst_up[i][isys][firstGood].Clone()
                h_bkg_tot_syst_dn = h_bkg_vecs_syst_dn[i][isys][firstGood].Clone()
                for icat in range(ncats):
                    if icat == firstGood: continue
                    if h_bkg_vecs_syst_up[i][isys][icat] != None:
                        h_bkg_tot_syst_up.Add(h_bkg_vecs_syst_up[i][isys][icat])
                    if h_bkg_vecs_syst_dn[i][isys][icat] != None:
                        h_bkg_tot_syst_dn.Add(h_bkg_vecs_syst_dn[i][isys][icat])
                if h_bkg_tot_syst_up.Integral() <= 0 or h_bkg_tot_syst_dn.Integral() <= 0:
                    print(h_bkg_tot_syst_up)
                h_bkg_tot_syst_up.Scale(h_bkg_tot.Integral()/h_bkg_tot_syst_up.Integral())
                h_bkg_tot_syst_dn.Scale(h_bkg_tot.Integral()/h_bkg_tot_syst_dn.Integral())
                h_bkg_tot_syst_up.Divide(h_bkg_tot)
                h_bkg_tot_syst_dn.Divide(h_bkg_tot)
                for ibin in range(1, h_bkg_tot.GetNbinsX()+1):
                    sysup = h_bkg_tot_syst_up.GetBinContent(ibin)-1
                    sysdn = h_bkg_tot_syst_dn.GetBinContent(ibin)-1
                    h_bkg_tot_syst_up.SetBinError(ibin, abs(sysup))
                    h_bkg_tot_syst_dn.SetBinError(ibin, abs(sysdn))
                    h_bkg_tot_syst_up.SetBinContent(ibin, 1)
                    h_bkg_tot_syst_dn.SetBinContent(ibin, 1)
                    # Temporary
                    syst_up[ibin-1] = ((syst_up[ibin-1])**2 + sysup**2)**0.5
                    syst_dn[ibin-1] = ((syst_dn[ibin-1])**2 + sysdn**2)**0.5

            systs[i] = syst_up # Temporary
            for ibin in range(len(syst_up)):
                systs[i][ibin] = max(syst_up[ibin], syst_dn[ibin])

    ## get data histograms
    if data==None:
        h_data = [None for i in plots]
    else:
        data_file = os.path.join(rootdir, data+".root")
        fid = ROOT.TFile(data_file)
        for i, pl in enumerate(plots):
            if i in skipList: continue
            vn = pl[0]
            if suffix != None: vn += suffix
            vni = vn+multisuf[0] if type(multisuf)==list else vn
            h = fid.Get("{0}/h_{1}".format(dirnames[0],vni))
            h_data.append( h )
            if type(h_data[-1])==type(ROOT.TObject()):
                raise Exception("No {0}/h_{1} histogram for {2}!".format(dirnames[0], vni, data))
            h_data[-1].SetDirectory(0)
            # handle the case with more than one directory
            for idir in range(1, len(dirnames)):
                vni = vn+multisuf[idir] if type(multisuf)==list else vn
                h_data[-1].Add(fid.Get("{0}/h_{1}".format(dirnames[idir],vni)))
        fid.Close()

    ## get signal histograms
    h_sig_vec = [[] for i in plots]
    sig_names = [[] for i in plots]
    if sig_points is not None:
        sig_file = os.path.join(rootdir, signame+".root")
        fid = ROOT.TFile(sig_file)
        for i, pl in enumerate(plots):
            if i in skipList: continue
            for spt in sig_points:
                vn = pl[0]+suffix+spt
                h = fid.Get("{0}/h_{1}".format(dirnames[0],vn))
                h_sig_vec[i].append( h )
                sig_names[i].append( utils.GetSampleName(signame+spt))
                if len(plots[i]) >= 5:
                    h_sig_vec[i][-1].Rebin(plots[i][4])
                if type(h_sig_vec[i][-1])==type(ROOT.TObject()):
                    raise Exception("No {0}/h_{1} histogram for {2}!".format(dirname, vn, signame))
                h_sig_vec[i][-1].SetDirectory(0)

                # handle the case with more than one directory
                for idir in range(1, len(dirnames)):
                    vni = vn[idir] if type(multisuf)==list else vn
                    h = fid.Get("{0}/h_{1}".format(dirnames[idir],vni))
                    h_sig_vec[i][-1].Add( h )
        fid.Close()

    # make the output directory if it doesn't exist
    if not os.path.isdir(os.path.join(output_dir,dirname+tag)):
        os.makedirs(os.path.join(output_dir,dirname+tag))

    # make all of the plots
    for i in range(len(plots)):
        if i in skipList: continue
        vn = plots[i][0]
        if suffix != None: vn += suffix
        userMin,userMax = None,None
        if plots[i][3]!=None:
            userMin = plots[i][3][0]
            userMax = plots[i][3][1]
        if len(plots[i]) >= 5:
            h_data[i].Rebin(plots[i][4])
        doOverflow = True
        if len(plots[i]) >= 6:
            doOverflow = plots[i][5]
        legcoord = (0.54,0.67,0.87,0.895)
        if len(plots[i]) >= 7:
            legpos = plots[i][6]
            if legpos.lower() == 'tl':
                legcoord = (0.14, 0.68, 0.5, 0.89)
            elif legcoord.lower() == 'bl':
                legcoord = (0.14, 0.18, 0.5, 0.39)
            elif legcoord.lower() == 'tm':
                legcoord = (0.34, 0.68, 0.7, 0.89)
            elif legcoord.lower() == 'bm':
                legcoord = (0.34, 0.18, 0.7, 0.39)
            elif legcoord.lower() == 'bl':
                legcoord = (0.14, 0.18, 0.5, 0.39)
            elif legcoord.lower() == 'br':
                legcoord = (0.58, 0.18, 0.92, 0.39)
            else:               # 'tr'
                legcoord = (0.54, 0.67, 0.87, 0.895)
            
        markerSize=0.8
        title = utils.GetCRName(dirname+tag, output_dir)
        xAxisTitle = h_bkg_vecs[i][0].GetXaxis().GetTitle()
        xAxisTitle = xAxisTitle.replace('E_{T}^{miss}', 'p_{T}^{miss}')
        unit = None
        if xAxisTitle == "":
            xAxisTitle = utils.GetVarName(vn)
            unit = utils.GetUnit(vn)
        elif '[GeV]' in xAxisTitle:
            unit = 'GeV'
            xAxisTitle = xAxisTitle.replace('[GeV]', '')

        subtitles = utils.GetSubtitles(dirname, vn+tag)
        if len(plots[i]) >= 7 and plots[i][6].lower() == 'tl':
            title = ''
            subtitles = []

        if h_data[i]!=None:
            if not scaleMC:
                subLegText = ["# Observed events: {ndata}"]
            elif type(scaleMC) == float:
                subLegText = ["MC scaled by {}".format(scaleMC),"# Data events: {ndata}"]
            else:
                subLegText = ["MC scaled by {datamcsf}","# Data events: {ndata}"]
        else:
            subLegText = None
        # subLegText = None
        if dogencat:
            sns = [utils.GetSampleName(s) for s in gencats]
            colors = [utils.GetColor(s) for s in gencats]
        else:
            sns = [utils.GetSampleName(s) for s in samples]
            colors = [utils.GetColor(s) for s in samples]
       
        pltname = vn if type(multisuf)!=list else vn + tag
        for ext in exts:
            # saveAs = os.path.join(output_dir,dirname+tag,"{0}_{1}.{2}".format(dirname,vn,ext))
            saveAs = os.path.join(output_dir,dirname+tag,"{1}.{2}".format(dirname,pltname,ext))
            ppm.plotDataMC(h_bkg_vecs[i], sns, h_data[i], doPause=False, xAxisTitle=xAxisTitle, lumi=lumi, lumiUnit='fb',
                           title=title, subtitles=subtitles, dataTitle=datatitle, xRangeUser=plots[i][2], isLog=plots[i][1], saveAs=saveAs,
                           scaleMCtoData=scaleMC, xAxisUnit=unit, userMin=userMin, userMax=userMax, doSort=False, customColors=colors,
                           markerSize=markerSize, titleSize=0.035, subtitleSize=0.033, legCoords=legcoord, legNCol=2,
                           subLegText=subLegText, subLegTextSize=0.036, doBkgError=True, doOverflow=doOverflow, cmsTextSize=0.04,
                           convertToPoisson=False, drawZeros=False, drawSystematicBand=drawSystematicBand, systematics=systs[i],
                           h_sig_vec=h_sig_vec[i], sig_names=sig_names[i], ratioType=ratioType, ratioTitle='Obs./Sim.',
                           yRangeUserRatio=ratioRange, doComposition=doComposition)

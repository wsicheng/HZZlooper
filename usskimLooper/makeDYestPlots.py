#!/usr/bin/env python
from __future__ import print_function, division

import os, sys
from math import sqrt
import ROOT as r

r.gROOT.SetBatch(1)
r.gStyle.SetOptStat('')

def moveOverFlowToLastBin1D(hist, newXmax=-1.0):
    nbin = hist.GetNbinsX();
    lastbin = nbin
    if newXmax > 0:
        lastbin = hist.GetXaxis().FindBin(newXmax) - 1

    if newXmax > 0 or hist.GetBinContent(nbin+1) > 0:
        err = r.Double();
        hist.SetBinContent(lastbin, hist.IntegralAndError(lastbin, -1, err));
        hist.SetBinError(lastbin, err);
        for ibin in range(lastbin+1, nbin+2):
            hist.SetBinContent(ibin, 0);
            hist.SetBinError(ibin, 0);

    if newXmax > 0:
        hist.GetXaxis().SetRangeUser(hist.GetXaxis().GetXmin(), newXmax)

    return hist


def DrawHeaderText(canvas, lumi=137.2, state="Preliminary"):
    tCMS = r.TLatex()
    tCMS.SetNDC(1)
    tCMS.SetTextFont(61)
    tCMS.SetTextAlign(11)
    tCMS.SetTextSize(0.0625)
    canvas.cd()
    tCMS.DrawLatex(canvas.GetLeftMargin(), 1.0-canvas.GetTopMargin()+0.012, "CMS")

    tplm = r.TLatex()
    tplm.SetNDC(1)
    tplm.SetTextFont(52)
    tplm.SetTextAlign(11)
    tplm.SetTextSize(0.048)
    canvas.cd()
    tplm.DrawLatex(canvas.GetLeftMargin()+0.110, 1.0-canvas.GetTopMargin()+0.012, state)

    if lumi:
        ttext = r.TLatex()
        ttext.SetNDC(1)
        ttext.SetTextFont(42)
        ttext.SetTextAlign(31)
        ttext.SetTextSize(0.048)
        canvas.cd()
        text = "{0} {1}^{{-1}} ({2} TeV)".format(lumi,"fb",13)
        ttext.DrawLatex(1.0-canvas.GetRightMargin()-0.01, 1.0-canvas.GetTopMargin()+0.012, text)


def countXrange(h, xmin, xmax):
    bmin = 0
    bmax = -1
    for i in range(1, h.GetNbinsX()+1):
        edge = h.GetBinLowEdge(i)
        if bmin == 0 and edge >= xmin:
            bmin = i
        if xmax > xmin and edge >= xmax:
            bmax = i-1
            break

    err = r.Double()
    res = h.IntegralAndError(bmin, bmax, err)
    return res, float(err)


# def getPoissonCountingConfidenceInterval_Frequentist(sw_total, swsq_total, CL):
#     quant = (1. - CL) / 2.0;
#     count = (swsq_total<=0. ? (sw_total==0. ? 0. : sw_total) : std::pow(sw_total, 2)/swsq_total);

#     vlow = (count == 0. ? 0. : r.Math.chisquared_quantile(quant, 2. * count) / 2.);
#     vhigh = r.Math.chisquared_quantile_c(quant, 2 * (count + 1.)) / 2.;
#     if count>0.:
#         vlow *= sw_total/count;
#         vhigh *= sw_total/count;

#     return vlow, vhigh

def makeNjetSF(year='2016', config=2):
    # indir = '../usskimLooper/output/v4_07_llgCR_{}'.format(year)
    indir = '../usskimLooper/output/v5_00_llgCR_{}'.format(year)

    # fdat = r.TFile('{}/{}.root'.format(indir, 'data_{}_llgskim'.format(year)))
    fdat = r.TFile('{}/{}.root'.format(indir, 'data_llgskim'))
    if config == 1:
        fbkg = r.TFile('{}/{}.root'.format(indir, 'allBkgs_llg'))
    else:
        fzg = r.TFile('{}/{}.root'.format(indir, 'ZGTo2LG'))
        fdy = r.TFile('{}/{}.root'.format(indir, 'DY'))
    
    hdat = fdat.Get('OffShell/h_njets_fullMET')
    if config == 1:
        hbkg = fbkg.Get('OffShell/h_njets_fullMET')
    elif config == 2:
        hdy  = fdy.Get('OffShell/h_njets_fullMET')
        hdat.Add(hdy, -1)
        hbkg = fzg.Get('OffShell/h_njets_fullMET')
    else:
        print('Unknown config')

    err = r.Double()
    ibc = hdat.IntegralAndError(3, -1, err)
    hdat.SetBinContent(3, ibc)
    hdat.SetBinError(3, err)

    ibc = hbkg.IntegralAndError(3, -1, err)
    hbkg.SetBinContent(3, ibc)
    hbkg.SetBinError(3, err)
    hdat.Divide(hbkg)


    print('| {} | ver{} '.format(year, config), end='')
    for ib in range(1, 4):
        print ('| {0:.3f} +- {1:.4f} '.format(hdat.GetBinContent(ib), hdat.GetBinError(ib)), end='')

    print('|')

def makeCombinedHists():

    year = '2018'
    indir = '../usskimLooper/output/v4_06_llgCR_{}'.format(year)
    flists = ['data_{}_llgskim'.format(year), 'ZGTo2LG', 'DY', 'others', 'allBkgs_llg']

    tarsuf = '_ll'
    srcsufs = ['_ee', '_mumu']

    # hlist = ['mtZZ', 'mtZZ_b1', 'njets', 'boson_pt', 'met', 'DjjVBF', 'min_dphijemt', 'dphi_boson_met', 'dphi_lljets_met']

    for fn in flists:
        f = r.TFile('{}/{}.root'.format(indir, fn), "UPDATE")
        dns = f.GetListOfKeys()
        for sdn in dns:
            sdn = sdn.GetName()
            ssuf0 = srcsufs[0]
            if not ssuf0 in sdn: continue
            tdn = sdn.replace(ssuf0, tarsuf)
            if tdn not in dns: tdir = f.mkdir(tdn)
            else: tdir = f.Get(tdn)
            sdir1 = f.Get('{}'.format(sdn))
            if not sdir1: print("ERROR: cannot find {} in file ".format(sdn))
            hlist = sdir1.GetListOfKeys()
            for hkey in hlist:
                hn = hkey.GetName()
                thn = hn.replace(ssuf0, tarsuf)
                htst = f.Get('{}/{}'.format(tdn, thn))
                if htst: continue
                hnew = f.Get('{}/{}'.format(sdn, hn)).Clone(thn)
                for ssuf in srcsufs[1:]:
                    hsi = f.Get('{}/{}'.format(sdn.replace(ssuf0, ssuf), hn.replace(ssuf0, ssuf)))
                    if not hsi:
                        if not 'h_el' in hn: print("ERROR: cannot find {}/{} in file ".format(sdn, hn).replace(ssuf0,ssuf))
                        continue
                    hnew.Add(hsi) 
                tdir.cd()
                hnew.Write()


def makeSubtraction():

    year = '2018'
    indir = '../usskimLooper/output/'
    skimver = 'v4_04'
    # phsuf = '_ee2j_flateta_rwgtd'
    phsuf = '2_ee2j_flateta_rwgtd_metge125'

    fdat = r.TFile('{0}{2}_phCR{3}_{1}/data_{1}_phskim.root'.format(indir, year, skimver, phsuf))
    fbkg = r.TFile('{0}{2}_phCR{3}_{1}/subtractor.root'.format(indir, year, skimver, phsuf))
    fout = r.TFile('{0}{2}_phCR{3}_{1}/DYestFromCR.root'.format(indir, year, skimver, phsuf), 'RECREATE')


    year = 'run2'
    indir = '../usskimLooper/output/fthists_'
    skimver = 'v4_06_dphi0p5_met140'
    fdat = r.TFile('{0}{2}_{1}/data.root'.format(indir, year, skimver, phsuf))
    fbkg = r.TFile('{0}{2}_{1}/subtractor.root'.format(indir, year, skimver, phsuf))
    fout = r.TFile('{0}{2}_{1}/DYestFromCR.root'.format(indir, year, skimver, phsuf), 'RECREATE')

    nfast = 'data single-#gamma'
    nfull = ''

    # hlist = ['mtZZ', 'mZZ', 'mtZZ_b1', 'mtZZ_b3', 'njets', 'boson_pt', 'met', 'DjjVBF', 'usmetrat' ]
    # hlist = ['mtZZ', 'mZZ', 'mtZZ_b1', 'mtZZ_b3', 'njets', 'boson_pt', 'met', 'DjjVBF', 'min_dphijemt', 'dphi_boson_met', 'dphi_lljets_met']
    hlist = ['mtZZ', 'mtZZ_b1', 'njets', 'boson_pt', 'met', 'DjjVBF', 'min_dphijemt', 'dphi_boson_met', 'dphi_lljets_met']

    for lsuf in ['_ll',]:
    # for metsuf in ['_fullMET', '_metge125',  '_metlt125']:
        for jsuf in ['_eq0j', '_eq1j', '_ge2j']:
            dn = 'OffShell'+jsuf+lsuf
            hdir = fdat.Get('{}'.format(dn))
            if not hdir:
                print("ERROR: cannot find {} in data ".format(dn))
            hlst = hdir.GetListOfKeys()

            for hkey in hlst:
                hn = hkey.GetName()
                inlist = False
                for hname in hlist:
                    if hname in hn: inlist = True

                if not inlist: continue

                hdat = fdat.Get('{}/{}'.format(dn, hn))
                hbkg = fbkg.Get('{}/{}'.format(dn, hn))
                if not hdat or not hbkg:
                    print("ERROR: cannot find {}/{} in either data or the subtractor".format(dn,hn))
                    continue

                dout = fout.Get(dn)
                if not dout: dout = fout.mkdir(dn)
                dout.cd()
                hout = hdat.Clone(hn)
                hsub = hbkg.Clone(hn+"_sub")
                hsub.Scale(-1.)
                hout.Add(hsub)
                for ibin in range(1, hout.GetNbinsX()):
                    ibc = hout.GetBinContent(ibin)
                    if ibc < 0.:
                        hout.SetBinContent(ibin, 0)
                        hout.SetBinError(ibin, 0)
                        obc = hdat.GetBinContent(ibin)
                        if obc > 0. and abs(ibc)/obc > 0.1 and 'mtZZ_b' in hn:
                            print( '[Subtraction] Setting output hist {} bin {} ({} to {}) from {} to 0, orig yield is {}'.format(hn, ibin, hout.GetBinLowEdge(ibin), hout.GetBinLowEdge(ibin+1), ibc, obc) )
                hout.Write()


def plotDYComparePlots(hdir, hname='met', samp='', rebin=1, y2=None, newrange=None, norm=True, **kwargs):

    sampsuf = '' if samp == '' else '_'+samp
    doscale = norm

    year  = kwargs.get('year', '2018')
    indir = '../usskimLooper/output/'
    phsuf = kwargs.get('phsuf', '_ee2j_rwgtd')
    llsuf = kwargs.get('llsuf', phsuf)
    # phsuf = '_ee2j_flateta_rwgtd'
    sdir = '{}/data{}_{}'.format(kwargs.get('sdir', 'dyest_v4_01_nosub'), phsuf, year)

    # y2 = True
    ffast = r.TFile('{0}{2}_phCR{3}_{1}/data_{1}_phskim.root'.format(indir, year, kwargs.get('skimver', 'v4_02'), phsuf))
    ffull = r.TFile('{0}{2}_2lSR{3}_{1}/data_{1}_llskim.root'.format(indir, year, kwargs.get('skimver', 'v4_02'), llsuf))
    if y2: fy2 = r.TFile('{0}{2}_phCR{2}_{1}/subtractor.root'.format(indir, year, kwargs.get('skimver', 'v4_02'), phsuf))
    nfast = 'data single-#gamma'
    nfull = 'data ee/#mu#mu'

    # gsuf = 'll' if 'rwgtd' in phsuf else 'gamma'
    # ffast = r.TFile('{}v4_01_metcrtd_phCR{}_{}/GJets.root'.format(indir, phsuf, year))
    # ffull = r.TFile('{}v4_01_metcrtd_2lSR_{}/DYJetsToLL_M-50.root'.format(indir, year))
    # nfast = '#gamma+Jets'
    # nfull = 'DY'
    # sdir = sdir.replace('data','closure')

    hfull = ffull.Get('{0}{2}_mumu/h_{1}{2}_mumu'.format(hdir,hname,sampsuf))
    hee = ffull.Get('{0}{2}_ee/h_{1}{2}_ee'.format(hdir,hname,sampsuf))
    lsuf = kwargs.get('lsuf', 'll')
    if lsuf == 'mumu':
        pass
    elif lsuf == 'ee':
        hfull = hee
    else: # 'll'
        hfull.Add(hee)

    gsuf = lsuf if 'rwgtd' in phsuf else 'gamma'
    hfast = ffast.Get('{0}{2}_{3}/h_{1}{2}_{3}'.format(hdir,hname,sampsuf,gsuf))
    if y2: hy2 = fy2.Get('{0}{2}_{3}/h_{1}{2}_{3}'.format(hdir,hname,sampsuf,gsuf))

    if not hfull: print( 'Cannot find {0}{2}/h_{1}{2} in {3}'.format(hdir,hname,sampsuf, 'ffull'))
    if not hfast: print( 'Cannot find {0}{2}/h_{1}{2} in {3}'.format(hdir,hname,sampsuf, 'ffast'))
    if y2 and not hy2: print( 'Cannot find {0}{2}/h_{1}{2} in {3}'.format(hdir,hname,sampsuf, 'ffast'))

    if y2:
        nhy2 = hy2.Clone(hy2.GetName()+"_clone")
        nhy2.Scale(-1.)
        hfast.Add(nhy2)
        for ibin in range(1, hfast.GetNbinsX()):
            ibc = hfast.GetBinContent(ibin)
            if ibc < 0.:
                print( '[{}] Setting bin {} ({} to {}) from {} to 0'.format(hname+sampsuf, ibin, hfast.GetBinLowEdge(ibin), hfast.GetBinLowEdge(ibin+1), ibc) )
                hfast.SetBinContent(ibin, 0)
                hfast.SetBinError(ibin, 0)

    rtype = '2#font[12]{l}/1#gamma'

    setrange = True if type(newrange) == list and len(newrange) > 1 else False
    newmax = newrange[1] if setrange else -1
    moveOverFlowToLastBin1D(hfull,newmax)
    moveOverFlowToLastBin1D(hfast,newmax)
    if y2: moveOverFlowToLastBin1D(hy2,newmax)

    if rebin != 1:
        hfull.Rebin(rebin)
        hfast.Rebin(rebin)

    if setrange:
        hfull.GetXaxis().SetRangeUser(newrange[0], newrange[1])
        hfast.GetXaxis().SetRangeUser(newrange[0], newrange[1])

    c0 = r.TCanvas('c0', 'c0', 800, 800)

    mainPad = r.TPad('1', '1', 0.0, 0.20, 1.0, 0.99)
    ratioPad = r.TPad('2', '2', 0.0, 0.02, 1.0, 0.24)

    r.SetOwnership(c0, False)
    r.SetOwnership(mainPad, False)
    r.SetOwnership(ratioPad, False)

    mainPad.SetTopMargin(0.08)
    mainPad.SetLeftMargin(0.12)
    mainPad.SetRightMargin(0.05)
    mainPad.SetBottomMargin(0.15)
    mainPad.SetLogy()
    mainPad.Draw()
    ratioPad.SetTopMargin(0.05)
    ratioPad.SetLeftMargin(0.12)
    ratioPad.SetRightMargin(0.05)
    ratioPad.SetBottomMargin(0.10)
    ratioPad.Draw()

    mainPad.cd()

    hfull.SetTitle('')
    hfull.GetYaxis().SetTitle('Events')
    hfull.GetYaxis().SetTitleOffset(0.75)
    hfull.GetYaxis().SetTitleSize(0.058)
    hfull.GetYaxis().SetTickSize(0.015)
    hfull.GetXaxis().SetTitleSize(0.04)

    hfull.SetLineWidth(2)
    hfast.SetLineWidth(2)

    hfull.SetMarkerSize(1.3)
    hfast.SetMarkerSize(1.3)

    hfull.SetMarkerStyle(25)
    hfast.SetMarkerStyle(26)

    hfull.SetLineColor(r.kGreen+3)
    hfast.SetLineColor(r.kRed)

    scale = hfull.Integral()/hfast.Integral()
    hncp = hfast.Clone("hfast_copy")
    if doscale:
        hfast.Scale(hfull.Integral()/hfast.Integral())

    hfull.Draw()
    hfast.Draw("same")


    leg = r.TLegend(0.58, 0.68, 0.92, 0.89)
    txtpos = (0.58, 0.64)
    selpos = (0.86, 0.56)

    if 'dphi' in hname:
        leg = r.TLegend(0.14, 0.18, 0.5, 0.39)
        txtpos = (0.14, 0.4)
        selpos = (0.16, 0.84)
    elif 'Efrac' in hname or 'eta' in hname:
        leg = r.TLegend(0.34, 0.18, 0.7, 0.39)
        txtpos = (0.34, 0.4)
    elif 'n' == hname[0] or 'Djj' in hname:
        leg = r.TLegend(0.14, 0.18, 0.5, 0.39)


    leg.AddEntry(hfast, nfast)
    leg.AddEntry(hfull, nfull)
    leg.Draw()

    selstr = ''
    if samp == 'eq0j':    selstr = '0j'
    elif samp == 'eq1j':  selstr = '1j'
    elif samp == 'eq2j':  selstr = '2j'
    elif samp == 'ge2j':  selstr = '#geq 2j'
    elif samp == 'geq1j': selstr = '#geq 1j'
    elif samp == 'geq2j': selstr = '#geq 2j'
    elif samp == 'vbf':   selstr = 'VBF sel.'

    atxt = r.TLatex();
    atxt.SetNDC();
    atxt.SetTextSize(0.04);
    # atxt.SetTextAlign(22);
    if doscale:
        atxt.DrawLatex(txtpos[0], txtpos[1], "scale ({0}) = {1:.3f}".format(rtype,scale));

    atxt.SetTextSize(0.08);
    atxt.DrawLatex(selpos[0], selpos[1], selstr);

    # DrawHeaderText(mainPad, 35.9 if '16' in samp else 41.5 if '17' in samp else 59.7 if '18' in samp else 137)
    # DrawHeaderText(mainPad, None, 'Simulation')
    DrawHeaderText(mainPad, None, 'Preliminary')

    ratioPad.cd()

    h_axis_ratio = r.TH1D('ratio_axis','', hfull.GetNbinsX(), hfull.GetXaxis().GetXmin(), hfull.GetXaxis().GetXmax())
    h_axis_ratio.GetYaxis().SetNdivisions(4)
    h_axis_ratio.GetYaxis().SetLabelSize(0.15)
    h_axis_ratio.GetYaxis().SetRangeUser(0, 2)
    # h_axis_ratio.GetYaxis().SetTitle(rtype)
    h_axis_ratio.GetYaxis().SetTitle('Ratio   ')
    # h_axis_ratio.GetYaxis().SetTitle('cms4 / OS')
    h_axis_ratio.GetYaxis().SetTitleOffset(0.25)
    h_axis_ratio.GetYaxis().SetTitleSize(0.18)
    h_axis_ratio.GetYaxis().SetTickLength(0.01)
    h_axis_ratio.GetXaxis().SetTickLength(0.07)
    h_axis_ratio.GetXaxis().SetTitleSize(0.0)
    h_axis_ratio.GetXaxis().SetLabelSize(0.0)
    if setrange: h_axis_ratio.GetXaxis().SetRangeUser(newrange[0], newrange[1])
    h_axis_ratio.Draw('axis')

    line = r.TF1('line1', '1', hfull.GetXaxis().GetXmin(), hfull.GetXaxis().GetXmax())
    line.SetLineStyle(2)
    line.SetLineColor(r.kGray+2)
    line.Draw('same')

    rat = hfast.Clone('ratio_gjetsvsdy')
    rat.Divide(hfull)
    rat.Draw("same")

    drawTFHist = False
    printweights = 'boson_pt_b' in hname and ('metlt125' in hname or 'metlt80' in hname)
    # printweights = printweights or ('boson_eta' in hname and 'metlt125' in hname and 'eq2j' in sampsuf)
    if 'flateta_rwgtd' in phsuf:
        printweights = False
    elif 'rwgtd_raweta' in phsuf:
        printweights = 'boson_aeta' in hname and 'metlt125' in hname and '2j' in sampsuf
    elif 'rwgtd' in phsuf:
        printweights = 'boson_aeta' in hname and 'metlt125' in hname and '2j' in sampsuf

    if 'eq0j' in sampsuf and 'boson_pt_b1' in hname: printweights = False
    if 'eq0j' not in sampsuf and 'boson_pt_b0' in hname: printweights = False

    if printweights:
        drawTFHist = True
        htf = hfull.Clone("h{}_{}_{}{}".format(year, hdir, hname, sampsuf))
        htf = hfull.Clone('ratio_2lvs1gamma')
        htf.Divide(hncp)
        if False:
            for i in range(1, htf.GetNbinsX()+1):
                print( "{0:.3e},".format(htf.GetBinLowEdge(i)), end='' )
            print( "" )
        jsuf = '_ee2j' if ('_ee2j' in phsuf or '_all2j' in phsuf) and '2j' in sampsuf  else sampsuf
        print( "vector<float> {} = {{".format(hname.replace('boson_','sf_V')+jsuf+'_data'+year[2:]), end='' )
        for i in range(1, htf.GetNbinsX()+1):
            ra = htf.GetBinContent(i)
            re = htf.GetBinError(i)
            if ra <= 0:
                htf.SetBinContent(i, 1.0)
                htf.SetBinError(i, 0.0)
            print( "{0:.3e},".format(ra if ra > 0 else 1.0), end='')
        print( "};" )

        # Print the error on transfer factors
        print( "vector<float> {} = {{".format(hname.replace('boson_','sferr_V')+jsuf+'_data'+year[2:]), end='' )
        for i in range(1, htf.GetNbinsX()+1):
            print( "{0:.3e},".format(htf.GetBinError(i)), end='')
        print( "};" )

    snsuf = "" if rebin==1 else "_rebin{}".format(rebin)

    savename = 'plots/{}/dyest_{}_{}{}{}.pdf'.format(sdir, hdir, hname, sampsuf, snsuf)
    os.system('mkdir -p plots/{}'.format(sdir))
    c0.SaveAs(savename)
    # print( 'Plot generated as ', savename)
    c0.SaveAs(savename.replace('.pdf','.png'))

    if drawTFHist:
        c1 = r.TCanvas('c1', 'c1', 800, 500)
        htf.SetMarkerStyle(1)
        htf.SetMarkerColor(r.kBlue)
        htf.SetLineColor(r.kBlue)
        if 'eta' in savename:
            htf.SetMaximum(2.0)
            htf.SetMinimum(0.0)
        else:
            htf.SetMaximum(0.2)
            htf.SetMinimum(0.0)
        htf.Draw()
        os.system('mkdir -p plots/{}'.format(sdir.replace('dyest','tfvals')))
        c1.SaveAs(savename.replace('dyest','tfvals'))
        c1.SaveAs(savename.replace('dyest','tfvals').replace('.pdf','.png'))
        frat = r.TFile('ratios.root', 'update')
        htf.Write()


def makeDYtestPlots():

    args = {
        'year'    : '2018',
        'sdir'    : 'dyest_v5_00_May9',
        # 'phsuf'   : '_all2jsel_rwgtd_raweta',
        'phsuf'   : '_all2jsel',
        'llsuf'   : '_all2jsel',
        'skimver' : 'v4_06',
        'lsuf'    : 'mumu',
    }

    # for metsuf in ['_metlt80', '_met80to125', '_metlt125']:
    for year in ['2016', '2017', '2018']:
        args['year'] = year
        for metsuf in [ '_metlt125']:
            for jsuf in ['eq0j', 'eq1j', 'ge2j']:
                plotDYComparePlots('OffShell', 'nvtxs'+metsuf, jsuf, **args)
                plotDYComparePlots('OffShell', 'boson_pt'+metsuf, jsuf, rebin=2, **args)
                plotDYComparePlots('OffShell', 'boson_eta'+metsuf, jsuf, rebin=1, **args)
                plotDYComparePlots('OffShell', 'boson_aeta'+metsuf, jsuf, rebin=1, **args)
                plotDYComparePlots('OffShell', 'boson_pt_b0'+metsuf, jsuf, **args)
                plotDYComparePlots('OffShell', 'boson_pt_b1'+metsuf, jsuf, **args)
                plotDYComparePlots('OffShell', 'boson_phi'+metsuf, jsuf, rebin=4, **args)
                plotDYComparePlots('OffShell', 'boson_mass'+metsuf, jsuf, **args)
                # plotDYComparePlots('OffShell', 'boson_phi_endcap'+metsuf, jsuf, rebin=4, **args)
                plotDYComparePlots('OffShell', 'mtZZ'+metsuf, jsuf, newrange=[150,600], rebin=4, **args)
                plotDYComparePlots('OffShell', 'mZZ'+metsuf,  jsuf, newrange=[150,600], rebin=4, **args)
                plotDYComparePlots('OffShell', 'met'+metsuf, jsuf, newrange=[0,450], rebin=2, **args)
                plotDYComparePlots('OffShell', 'dphi_boson_met'+metsuf, jsuf, **args)
                plotDYComparePlots('OffShell', 'dphi_lljets_met'+metsuf, jsuf, **args)

            # for jsuf in ['ge2j']:
            #     plotDYComparePlots('OffShell', 'DjjVBF'+metsuf, jsuf, rebin=4, **args)


if __name__ == '__main__':

    makeDYtestPlots()
    # makeSubtraction()
    # makeCombinedHists()

    makeNjetSF('2016', 1)
    makeNjetSF('2016', 2)
    makeNjetSF('2017', 1)
    makeNjetSF('2017', 2)
    makeNjetSF('2018', 1)
    makeNjetSF('2018', 2)


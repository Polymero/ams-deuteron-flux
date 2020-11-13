# A SIMPLE CUT ANALYSIS (See PyROOT Geo-magnetic Cut-off.ipynb)
# Sebastiaan Venendaal
# Created on        26-10-20
# Last edited on    27-10-20

# Imports
import numpy as np
import ROOT

# ------------------------------------------------------------------------------
# DATA
# ------------------------------------------------------------------------------
# Output directory for images
out_path = './CutAnalysis/'
# Include header file
ROOT.gInterpreter.Declare('#include "./Ntp.h"')
# Make a Tree Chain
chain = ROOT.TChain('Simp')
# Fill Chain with all available Runs
chain.Add('../Simp.root')
# Print number of entries
chain_events = chain.GetEntries()
print('Data Entries in TChain:', chain_events)

# ------------------------------------------------------------------------------
# SELECTION CUTS
# ------------------------------------------------------------------------------
# Cut configuration
lay_crit = '>=' # criterion symbol for bad charge status in tracker layers
geo_crit = '1.2' # criterion strength of geo-magnetic cut-off (standard convention is 1.2)
# General
CGa = ROOT.TCut('tof_beta > 0') # down-going particles
CGb = ROOT.TCut('trk_rig > 0 && trk_rig <= 22') # deuteron identification range
CGc = ROOT.TCut('rich_beta > 0.95 && rich_beta < 0.98') # small beta slice (optional)
CG = CGa + CGb
# Quality
CQa = ROOT.TCut('trk_chisqn[0] < 10 && trk_chisqn[1] < 10 && trk_chisqn[0] > 0 && trk_chisqn[1] > 0') # well-reconstructed particles
CQb = ROOT.TCut('status % 10 == 1') # single particle events
CQ = CQa + CQb
# Charge
CZa = ROOT.TCut('trk_q_inn > 0.80 && trk_q_inn < 1.30') # single charge particles
CZb = ROOT.TCut('trk_q_lay[4]{0}0 && trk_q_lay[1]{0}0 && trk_q_lay[2]{0}0 && trk_q_lay[3]{0}0 &&'.format(lay_crit) +
                'trk_q_lay[5]{0}0 && trk_q_lay[6]{0}0 && trk_q_lay[7]{0}0 && trk_q_lay[8]{0}0 &&'.format(lay_crit) +
                'trk_q_lay[0]{0}0'.format(lay_crit)) # no bad charge status throughout tracker
CZ = CZa + CZb
# Geo-magnetic cut-off
CR = ROOT.TCut('trk_rig > {0}*cf'.format(geo_crit)) # geo-magnetic cut-off
# TOF - RICH Consistency
CC = ROOT.TCut('abs(tof_beta-rich_beta)/tof_beta < 0.05')

# ------------------------------------------------------------------------------
# PLOT FUNCTION
# ------------------------------------------------------------------------------

ROOT.gStyle.SetOptTitle(0)

def cut_plot(par, title='', xlabel='', ylabel='events', stats=False, bins=200, xrange=(0,0), yrange=(0,0), log_y=True, prop=False, lx=0.65, ly=0.6):
    # Create Canvas
    if title == '':
        title = par
    canvas = ROOT.TCanvas("{}".format(title),"{}".format(title), 1600, 1000)
    # Draw Histogram and store in ROOT.???
    xrange = str(xrange)[1:-1]
    chain.Draw('{} >> {}({}, {})'.format(par, title+'1', bins, xrange), '', '')
    chain.Draw('{} >> {}({}, {})'.format(par, title+'2', bins, xrange), CG, 'SAME')
    chain.Draw('{} >> {}({}, {})'.format(par, title+'3', bins, xrange), CG+CQ, 'SAME')
    chain.Draw('{} >> {}({}, {})'.format(par, title+'4', bins, xrange), CG+CQ+CZ, 'SAME')
    chain.Draw('{} >> {}({}, {})'.format(par, title+'5', bins, xrange), CG+CQ+CZ+CR, 'SAME')
    chain.Draw('{} >> {}({}, {})'.format(par, title+'6', bins, xrange), CG+CQ+CZ+CR+CC, 'SAME')
    # Get TH1F objects
    hist1 = ROOT.gROOT.FindObject(title+'1')
    hist2 = ROOT.gROOT.FindObject(title+'2')
    hist3 = ROOT.gROOT.FindObject(title+'3')
    hist4 = ROOT.gROOT.FindObject(title+'4')
    hist5 = ROOT.gROOT.FindObject(title+'5')
    hist6 = ROOT.gROOT.FindObject(title+'6')
    # Coloring
    hist1.SetLineColor(ROOT.kBlue); hist1.SetLineWidth(3)
    hist2.SetLineColor(ROOT.kRed); hist2.SetLineWidth(3)
    hist3.SetLineColor(ROOT.kGreen); hist3.SetLineWidth(3)
    hist4.SetLineColor(ROOT.kOrange); hist4.SetLineWidth(3)
    hist5.SetLineColor(ROOT.kBlack); hist5.SetLineWidth(5)
    hist6.SetLineColor(ROOT.kMagenta); hist6.SetLineWidth(3)
    # Labelling
    # hist1.SetTitle(title) # Title is bugging, so just get rid of it (won't use it in thesis anyway)
    hist1.GetXaxis().SetTitle(xlabel)
    hist1.GetYaxis().SetTitle(ylabel)
    hist1.SetStats(stats)
    # Axis
    if yrange != (0,0):
        hist1.SetAxisRange(yrange[0], yrange[1], "Y")
    canvas.SetLogy(log_y)
    # Legend
    legend = ROOT.TLegend(lx, ly,lx+.20, ly+.15)
    legend.AddEntry(hist1, 'None')
    legend.AddEntry(hist2, 'CG')
    legend.AddEntry(hist3, 'CG+CQ')
    legend.AddEntry(hist4, 'CG+CQ+CZ')
    legend.AddEntry(hist5, 'CG+CQ+CZ+CR')
    legend.AddEntry(hist6, 'CG+CQ+CZ+CR+CC')
    legend.SetLineWidth(0)
    #legend.SetTextSize(10) # This line makes text in TLegend disappear !!!
    legend.Draw() # TLegend can be drawn here but has to returned as well !!!
    # Proportion of events
    if prop is True:
        he1 = hist1.GetEntries(); print('None: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
        he1 = hist2.GetEntries(); print('CG: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
        he1 = hist3.GetEntries(); print('CG+CQ: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
        he1 = hist4.GetEntries(); print('CG+CQ+CZ: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
        he1 = hist5.GetEntries(); print('CG+CQ+CZ+CR: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
        he1 = hist6.GetEntries(); print('CG+CQ+CZ+CR+CC: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
    # Delete histograms
    # ROOT.hist1.SetDirectory(0)
    # ROOT.hist2.SetDirectory(0)
    # ROOT.hist3.SetDirectory(0)
    # ROOT.hist4.SetDirectory(0)
    # ROOT.hist5.SetDirectory(0)
    # Return canvas and legend elements
    return canvas, legend

# ------------------------------------------------------------------------------
# CUT PLOTS
# ------------------------------------------------------------------------------
# TOF Beta (prop = True)
c, l = cut_plot('tof_beta', title='TOF Beta', xlabel='beta (v/c)', xrange=(.01,2.), prop=True)
c.Print(out_path + 'TOF_Beta.png')
# RICH Beta
c, l = cut_plot('rich_beta', title='RICH Beta', xlabel='beta (v/c)', xrange=(.74,1.05), lx=.2)
c.Print(out_path + 'RICH_Beta.png')
# Tracker Rigidity
c, l = cut_plot('trk_rig', title='Tracker Rigidity', xlabel='R [GV]', xrange=(0,22), yrange=(1,1e6))
c.Print(out_path + 'Tracker_Rigidity.png')
# Inner Tracker Charge
c, l = cut_plot('trk_q_inn', title='Inner Tracker Charge', xlabel='Z [e]', xrange=(0.3,3.0), yrange=(1,1e5), lx=.7)
c.Print(out_path + 'Inner_Tracker_Charge.png')
# First Tracker Layer Charge
c, l = cut_plot('trk_q_lay[0]', title='First Tracker Layer Charge', xlabel='Z [e]', xrange=(0.3,3.0), yrange=(1,3e4), lx=.7)
c.Print(out_path + 'First_Tracker_Layer_Charge.png')
# Bending Plane Chi-squared
c, l = cut_plot('trk_chisqn[1]', title='Bending Plane Chi-squared', xlabel='chi-squared', xrange=(-.5,10.5), yrange=(1,3e5))
c.Print(out_path + 'Bending_Plane_Chisquared.png')
# TOF Mass
c, l = cut_plot('trk_q_inn * trk_rig * TMath::Sqrt(1/tof_beta/tof_beta - 1)', title='TOF Mass Log',
                xlabel='m [GeV/c^2]', xrange=(.1,2.5), yrange=(1, 3750), log_y=True, ly=.15)
c.Print(out_path + 'TOF_Mass_Log.png')
c, l = cut_plot('trk_q_inn * trk_rig * TMath::Sqrt(1/tof_beta/tof_beta - 1)', title='TOF Mass',
                xlabel='m [GeV/c^2]', xrange=(.1,2.5), yrange=(0, 4000), log_y=False)
c.Print(out_path + 'TOF_Mass.png')
# RICH Mass
c, l = cut_plot('trk_q_inn * trk_rig * TMath::Sqrt(1/rich_beta/rich_beta - 1)', title='RICH Mass Log',
                xlabel='m [GeV/c^2]', xrange=(.1,2.5), yrange=(1, 4250), log_y=True, ly=.7)
c.Print(out_path + 'RICH_Mass_Log.png')
c, l = cut_plot('trk_q_inn * trk_rig * TMath::Sqrt(1/rich_beta/rich_beta - 1)', title='RICH Mass',
                xlabel='m [GeV/c^2]', xrange=(.1,2.5), yrange=(0, 5000), log_y=False)
c.Print(out_path + 'RICH_Mass.png')
# Livetime
c, l = cut_plot('lf', title='Livetime', xlabel='t [s]', xrange=(0,1), log_y=True, lx=.2);
c.Print(out_path + 'Livetime.png')
# TOF - RICH Consistency
c, l = cut_plot('abs(tof_beta-rich_beta)/tof_beta', title='TOF - RICH Beta Consistency',
                xlabel='abs(tof_beta-rich_beta)/tof_beta', xrange=(0,0.3), yrange=(1, 1e4), log_y=True)
c.Print(out_path + 'TOF_RICH_Beta_Consistency.png')

# ------------------------------------------------------------------------------
# YET TO FIX: ALL BELOW
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# GEO-MAGNETIC CUT-OFF PLOTS
# ------------------------------------------------------------------------------
# 2D RGeo VS Rigidity
c, l =ROOT.TCanvas("beta_rig","beta_rig", 800, 500)
chain.Draw('{0}*cf:trk_rig >> hist(200)'.format(geo_crit), CG+CQ+CZ, 'COLZ')
# Labelling
ROOT.hist.GetXaxis().SetTitle('R [GV]')
ROOT.hist.GetYaxis().SetTitle('{0}*Rgeo [GV]'.format(geo_crit))
ROOT.hist.SetTitle('Rgeo vs R (CG,CQ,CZ)')
ROOT.hist.SetStats(0)
# Axis
c.SetLogy(False)
c.SetLogx(False)


# EFFECT OF DIFFERENT RGeo CRITERIA
bins = 200 # number of bins
rang = '0.01, 22' # x-axis plot range
# Canvas and Histrograms
c, l =ROOT.TCanvas("Q_lay1_inn","Q_lay1_inn", 800, 500)
chain.Draw('trk_rig >> hist1({}, {})'.format(bins, rang), CG+CQ+CZ, '')
chain.Draw('trk_rig >> hist2({}, {})'.format(bins, rang), CG+CQ+CZ+ROOT.TCut('trk_rig>0.8*cf'), 'SAME')
chain.Draw('trk_rig >> hist3({}, {})'.format(bins, rang), CG+CQ+CZ+ROOT.TCut('trk_rig>1.0*cf'), 'SAME')
chain.Draw('trk_rig >> hist4({}, {})'.format(bins, rang), CG+CQ+CZ+ROOT.TCut('trk_rig>1.2*cf'), 'SAME')
chain.Draw('trk_rig >> hist5({}, {})'.format(bins, rang), CG+CQ+CZ+ROOT.TCut('trk_rig>1.4*cf'), 'SAME')
# Coloring
ROOT.hist1.SetLineColor(ROOT.kBlue)
ROOT.hist2.SetLineColor(ROOT.kRed)
ROOT.hist3.SetLineColor(ROOT.kGreen)
ROOT.hist4.SetLineColor(ROOT.kOrange)
ROOT.hist5.SetLineColor(ROOT.kBlack)
ROOT.hist4.SetLineWidth(3)
# Labelling
ROOT.hist1.GetXaxis().SetTitle('R [GV]')
ROOT.hist1.GetYaxis().SetTitle('events')
ROOT.hist1.SetTitle('Rigidity effect Rgeo Cut (CG, CQ, CZ)')
ROOT.hist1.SetStats(False)
# Axis
ROOT.hist1.SetAxisRange(1,1e4,"Y")
c.SetLogy(True)
# Legend
legend = ROOT.TLegend(0.65 ,0.6 ,0.85 ,0.75)
legend.AddEntry(ROOT.hist1, 'No Rgeo')
legend.AddEntry(ROOT.hist2, '0.8 Rgeo')
legend.AddEntry(ROOT.hist3, '1.0 Rgeo')
legend.AddEntry(ROOT.hist4, '1.2 Rgeo')
legend.AddEntry(ROOT.hist5, '1.4 Rgeo')
legend.SetLineWidth(0)
legend.SetTextSize(10)
legend.Draw("same")

# Proportion of events
he1 = ROOT.hist1.GetEntries(); print('\nNo Rgeo: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
he1 = ROOT.hist2.GetEntries(); print('0.8 Rgeo: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
he1 = ROOT.hist3.GetEntries(); print('1.0 Rgeo: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
he1 = ROOT.hist4.GetEntries(); print('1.2 Rgeo: ', he1, '({:.1f} %)'.format(he1/chain_events*100))
he1 = ROOT.hist5.GetEntries(); print('1.4 Rgeo: ', he1, '({:.1f} %)'.format(he1/chain_events*100))

# ------------------------------------------------------------------------------
# OTHER PLOTS
# ------------------------------------------------------------------------------
# TOF Mass VS Rigidity
c, l =ROOT.TCanvas("mass_rig","mass_rig", 800, 500)
chain.Draw('trk_q_inn * trk_rig * TMath::Sqrt(1/tof_beta/tof_beta - 1):trk_rig >> hist(200)'.format(geo_crit), CG+CQ+CZ+CR, 'COLZ')
# Labelling
ROOT.hist.GetXaxis().SetTitle('R [GV]')
ROOT.hist.GetYaxis().SetTitle('m [GeV/c^2]')
ROOT.hist.SetTitle('TOF Mass vs Rigidity (CG,CQ,CZ,CR)')
ROOT.hist.SetStats(0)
# Axis
c.SetLogy(False)
c.SetLogx(False)


# RICH Mass VS Rigidity
c, l =ROOT.TCanvas("rmass_rig","rmass_rig", 800, 500)
chain.Draw('trk_q_inn * trk_rig * TMath::Sqrt(1/rich_beta/rich_beta - 1):trk_rig >> hist(200)'.format(geo_crit), CG+CQ+CZ+CR, 'COLZ')
# Labelling
ROOT.hist.GetXaxis().SetTitle('R [GV]')
ROOT.hist.GetYaxis().SetTitle('m [GeV/c^2]')
ROOT.hist.SetTitle('RICH Mass vs Rigidity (CG,CQ,CZ,CR)')
ROOT.hist.SetStats(0)
# Axis
c.SetLogy(False)
c.SetLogx(False)


# Chi-squared VS Rigidity
c, l =ROOT.TCanvas("chi_rig","chi_rig", 800, 500)
chain.Draw('trk_chisqn[1]:trk_rig >> hist(200)'.format(geo_crit), CG+CQ+CZ+CR, 'COLZ')
# Labelling
ROOT.hist.GetXaxis().SetTitle('R [GV]')
ROOT.hist.GetYaxis().SetTitle('chi-squared')
ROOT.hist.SetTitle('Chi-squared (y-plane) vs Rigidity (CG,CQ,CZ,CR)')
ROOT.hist.SetStats(0)
# Axis
c.SetLogy(False)
c.SetLogx(False)


# Chi-squared VS TOF beta
c, l =ROOT.TCanvas("chi_beta","chi_beta", 800, 500)
chain.Draw('trk_chisqn[1]:tof_beta >> hist(200)'.format(geo_crit), CG+CQ+CZ+CR, 'COLZ')
# Labelling
ROOT.hist.GetXaxis().SetTitle('beta (v/c)')
ROOT.hist.GetYaxis().SetTitle('chi-squared')
ROOT.hist.SetTitle('Chi-squared (y-plane) vs TOF Beta (CG,CQ,CZ,CR)')
ROOT.hist.SetStats(0)
# Axis
c.SetLogy(False)
c.SetLogx(False)

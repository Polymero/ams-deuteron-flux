# A SIMPLE CUT EFFICIENCY CALCULATOR (See PyROOT Cut Efficiency.ipynb)
# Sebastiaan Venendaal
# Created on        01-11-20
# Last edited on    01-11-20

# Imports
import pandas as pd
import numpy as np
import ROOT

# ------------------------------------------------------------------------------
# DATA
# ------------------------------------------------------------------------------
# Output directory for images
out_path = './CutEfficiency/'
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

# CUTS DataFrame
cut_df = pd.DataFrame(columns=['Name', 'Cut', 'Instrument', 'Description'])
# Append all single CUTS
cut_df = cut_df.append({'Name': 'Crig', 'Cut': ROOT.TCut('trk_rig > 0 && trk_rig <= 22'),
          'Instrument': 'General Cut', 'Description': 'Deuteron identification range'}, ignore_index=True)
cut_df = cut_df.append({'Name': 'Cpar', 'Cut': ROOT.TCut('status % 10 == 1'),
          'Instrument': 'Combi', 'Description': 'Single particle events'}, ignore_index=True)
cut_df = cut_df.append({'Name': 'Ccon', 'Cut': ROOT.TCut('abs(tof_beta-rich_beta)/tof_beta < 0.05'),
          'Instrument': 'TOF', 'Description': 'Consistancy between TOF and RICH beta measurement'}, ignore_index=True)
cut_df = cut_df.append({'Name': 'Cbet', 'Cut': ROOT.TCut('tof_beta > 0'),
          'Instrument': 'TOF', 'Description': 'Down-going particles'}, ignore_index=True)
cut_df = cut_df.append({'Name': 'Cchi', 'Cut': ROOT.TCut('trk_chisqn[0] < 10 && trk_chisqn[1] < 10 && trk_chisqn[0] > 0 && trk_chisqn[1] > 0'),
          'Instrument': 'Tracker', 'Description': 'Well-reconstructed particle tracks'}, ignore_index=True)
cut_df = cut_df.append({'Name': 'Cinn', 'Cut': ROOT.TCut('trk_q_inn > 0.80 && trk_q_inn < 1.30'),
          'Instrument': 'Tracker', 'Description': 'Single charge particles'}, ignore_index=True)
cut_df = cut_df.append({'Name': 'Clay', 'Cut': ROOT.TCut('trk_q_lay[4]{0}0 && trk_q_lay[1]{0}0 && trk_q_lay[2]{0}0 && trk_q_lay[3]{0}0 &&'.format(lay_crit) +
                                     'trk_q_lay[5]{0}0 && trk_q_lay[6]{0}0 && trk_q_lay[7]{0}0 && trk_q_lay[8]{0}0 &&'.format(lay_crit) +
                                     'trk_q_lay[0]{0}0'.format(lay_crit)),
          'Instrument': 'Tracker', 'Description': 'No bad charge status throughout tracker'}, ignore_index=True)
cut_df = cut_df.append({'Name': 'Cgeo', 'Cut': ROOT.TCut('trk_rig > {0}*cf'.format(geo_crit)),
          'Instrument': 'Tracker', 'Description': 'Geo-magnetic cut-off'}, ignore_index=True)
# Print CUTS for visual check
print(cut_df)

# ------------------------------------------------------------------------------
# PLOT FUNCTION
# ------------------------------------------------------------------------------

ROOT.gStyle.SetOptTitle(0)

def eff_plot(cut, title='', xlabel='Rigidity [GV]', ylabel='efficiency', bins=50, lx=0.65, ly=0.2):
    # Create Canvas
    if title == '':
        title = par
    canvas = ROOT.TCanvas("{}".format(title),"{}".format(title), 1600, 1000)
    # Current cut
    current_inst = cut_df[cut_df['Name'] == cut]['Instrument'].values[0]
    current_cut  = cut_df[cut_df['Name'] == cut]['Cut'].values[0]
    # Other 'independent' cuts
    base_cut = ROOT.TCut(''); effi_cut = ROOT.TCut('')
    base_df = cut_df[(cut_df['Instrument'] != current_inst)]
    for c in base_df['Cut'].values:
        base_cut += c
    effi_df = sub_df.append(cut_df[cut_df['Name'] == 'Cchi'])
    for c in effi_df['Cut'].values:
        effi_cut += c
    # Draw Histogram and store in ROOT.???
    chain.Draw('trk_rig >> {}({}, 0, 22)'.format(title+'1', bins), base_cut, '')
    chain.Draw('trk_rig >> {}({}, 0, 22)'.format(title+'2', bins), effi_cut, 'SAME')
    # Get TH1F objects
    hist1 = ROOT.gROOT.FindObject(title+'1')
    hist2 = ROOT.gROOT.FindObject(title+'2')
    # Create ratio histogram
    hist3 = hist1.CLone(title+'3')
    hist3.Divide(hist2)
    # Axis
    hist1.SetAxisRange(0, 1, "Y")
    # Coloring
    hist3.SetLineColor(ROOT.kGreen); hist3.SetLineWidth(3)
    # Labelling
    # hist1.SetTitle(title) # Title is bugging, so just get rid of it (won't use it in thesis anyway)
    hist3.GetXaxis().SetTitle(xlabel)
    hist3.GetYaxis().SetTitle(ylabel)
    hist3.SetStats(stats)
    # Legend
    legend = ROOT.TLegend(lx, ly,lx+.20, ly+.15)
    legend.AddEntry(hist3, title)
    legend.SetLineWidth(0)
    #legend.SetTextSize(10) # This line makes text in TLegend disappear !!!
    legend.Draw() # TLegend can be drawn here but has to returned as well !!!
    # Return canvas and legend elements
    return canvas, legend

# ------------------------------------------------------------------------------
# CUT EFFICIENCY PLOTS
# ------------------------------------------------------------------------------
# TOF Beta cut
c, l = eff_plot()
c.Print(out_path + 'Downward Particle Cut (TOF Beta).png')

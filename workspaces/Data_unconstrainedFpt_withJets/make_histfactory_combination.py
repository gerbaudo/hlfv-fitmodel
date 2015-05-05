#!/usr/bin/env python
import os
import sys
import argparse
import xml.dom.minidom
import ROOT

#if len(sys.argv)<3:
#    print "Usage: %s <mass> <lumi> <lumierr> <output> [input1 input2 ...]" % os.path.basename(sys.argv[0])
#    exit()

parser = argparse.ArgumentParser()
parser.add_argument('-binlo', default='0')
parser.add_argument('-binhi', default='1')
parser.add_argument('-prefix', default='ws')
parser.add_argument('-makechecks', default=False, action='store_true')
parser.add_argument('-combinelumi', action='store_true')
parser.add_argument('lumi')
parser.add_argument('lumierr')
parser.add_argument('output')
parser.add_argument('ifiles',nargs='+')
args = parser.parse_args()
sys.argv=[]

lumival=float(args.lumi)
lumierr=float(args.lumierr)
binlo=int(args.binlo)
binhi=int(args.binhi)

syst=[]
norm=[]
stat=[]
cons=[]

for f in args.ifiles:
    doc = xml.dom.minidom.parse(f)
    chan= doc.getElementsByTagName("Channel")[0].attributes["Name"].value
    for e in doc.getElementsByTagName("StatError"):
        if e.attributes["Activate"].value == "True": 
            for i in range(20):
                stat.append("gamma_stat_%s_bin_%i"%(chan,i))
    for e in doc.getElementsByTagName("OverallSys"):
        if e.attributes["Name"].value not in syst: syst.append(e.attributes["Name"].value)
    for e in doc.getElementsByTagName("HistoSys"):
        if e.attributes["Name"].value not in syst: syst.append(e.attributes["Name"].value)
    for e in doc.getElementsByTagName("NormFactor"):
        if e.attributes["Name"].value == "mu_BR_htm": continue
        if e.attributes["Const"].value == "True":
            if e.attributes["Name"].value not in cons: cons.append(e.attributes["Name"].value)
        else:
            if e.attributes["Name"].value not in norm: norm.append(e.attributes["Name"].value)
    for e in doc.getElementsByTagName("ShapeFactor"):
#        if e.attributes["Const"].value == "True":
#            if e.attributes["Name"].value not in cons: cons.append(e.attributes["Name"].value)
#        else:
            if e.attributes["Name"].value not in norm: norm.append(e.attributes["Name"].value)
        
if args.combinelumi == True:
    cons.reverse()
    cons.append("Lumi")
    cons.reverse()

print "************************* in make_combination  *************************************"
    
for i in range(len(syst)):
    if syst[i].lower()=="lumi": continue
    if not syst[i].startswith("alpha_"):
        syst[i]="alpha_"+syst[i]
        
if not os.path.exists("results"):
    os.mkdir("results")

print "Writing: %s" % args.output 
fw = open(args.output,"w+")
fw.write("<!DOCTYPE Combination SYSTEM 'HistFactorySchema.dtd'>\n")
smode=''
if ROOT.gROOT.GetVersionInt()<53400: smode='Mode="comb"'
fw.write('<Combination OutputFilePrefix="./results/%s" %s>\n' % (args.prefix,smode))

for f in args.ifiles:
    fw.write('  <Input>%s</Input>\n'%f)


print "Generate All Systematics measurement: scale factors are constant"
fw.write('  <Measurement Name="AllSYS" Lumi="%.4f" LumiRelErr="%.4f" ExportOnly="True" >\n' % (lumival,lumierr))
fw.write('    <POI>mu_BR_htm</POI>\n')
if len(norm)>0:
    fw.write('    <ParamSetting Const="True">%s</ParamSetting>\n' % ' '.join(cons))
fw.write('  </Measurement>\n')


#print "Generate No Systematics measurement: scale factors, systematics and statistical parameters are constant"
#fw.write('  <Measurement Name="NoSYS" Lumi="%.4f" LumiRelErr="%.4f" ExportOnly="True" >\n' % (lumival,lumierr))
#fw.write('    <POI>SigXsecOverSM</POI>\n')
#if len(norm+cons+syst+stat)>0:
#    if args.combinelumi == False:
#        fw.write('    <ParamSetting Const="True">%s</ParamSetting>\n' % ' '.join(cons+syst+stat))
#    else:
#        nsyst = []
#        for s in syst:
#            if "lumi" in s.lower(): continue
#            nsyst.append(s)
#        fw.write('    <ParamSetting Const="True">%s</ParamSetting>\n' % ' '.join(cons+nsyst+stat))
#fw.write('  </Measurement>\n')


#print "Generate No Statistics measurement: scale factors and statistical parameters are constant"
#fw.write('  <Measurement Name="NoStat" Lumi="%.4f" LumiRelErr="%.4f" ExportOnly="True" >\n' % (lumival,lumierr))
#fw.write('    <POI>SigXsecOverSM</POI>\n')
#if len(norm+stat)>0:
#    fw.write('    <ParamSetting Const="True">%s</ParamSetting>\n' % ' '.join(cons+stat))
#fw.write('  </Measurement>\n')


if args.makechecks:
    for i in range(len(syst)-1):
        fw.write('  <Measurement Name="Check_%02i" Lumi="%.4f" LumiRelErr="%.4f" ExportOnly="True" >\n' % (i,lumival,lumierr))
        fw.write('    <POI>mu_BR_htm</POI>\n')
        fw.write('    <!-- Only %s is floating -->\n' % syst[i])
        fw.write('    <ParamSetting Const="True">%s</ParamSetting>\n' % ' '.join(syst[0:i]+syst[i+1:]))
        fw.write('  </Measurement>\n')


fw.write('</Combination>\n')
fw.close()

#Carlos: you also want to add a PDF for the constraint terms 
#for s in sorted(syst):
#    fw.write('    <ConstraintTerm Type="LogNormal">%s</ConstraintTerm>\n'%s)

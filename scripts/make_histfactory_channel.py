#!/usr/bin/env python
#
# Make HistFactory channel file for LFV
#
# prod_i { P(obs_i, prod_j{(1-d_j*e_j)} * sum_j{ bg_ji } +
#		    prod_k{(1-d_k*e_k)} * sum_j{ mu*sg_ji} )}
#
# Carlos.Solans@cern.ch

import os
import sys
import ROOT
import argparse

def isNumeric(n):
	try:
		i = float(n)
	except ValueError, TypeError:
		return False
	return True

def checkExistance(fr,n1):
	fr = ROOT.TFile(fr,'READONLY')
	h1 = fr.Get(n1)
	if h1==None: return False
	return True

def dirname(s):
	s=os.path.dirname(s)
	if s==None: return ""
	elif len(s)>1 and s[0]=="/": return s[1:]
	return s

memfile={}
def checkIdentical(fname,hnom,hvhi,hvlo):
	sys.argv=[]
	if fname not in memfile:
		fr = ROOT.TFile(fname,'READONLY')
		memfile[fname] = fr
	else:
		fr = memfile[fname]
	h0 = fr.Get(hnom)
	h1 = fr.Get(hvhi)
	h2 = fr.Get(hvlo)
	print h0.GetNbinsX(),hnom,hvhi,hvlo
    #print hnom,hvhi,hvlo
	for i in xrange(h0.GetNbinsX()+2):
		if(h0.GetBinContent(i)-h1.GetBinContent(i)>1e-6): return False
		if(h0.GetBinContent(i)-h2.GetBinContent(i)>1e-6): return False
	return True

parser = argparse.ArgumentParser()
parser.add_argument('-sig', nargs='?', default=[], help='sig histos', metavar='name:path', action='append')
parser.add_argument('-wrsig', nargs='?', default=[], help='wrong sig histos', metavar='name:path', action='append')
parser.add_argument('-bkg', nargs='?', default=[], help='bkg histos', metavar='name:path', action='append')
parser.add_argument('-obs', help='observed', metavar='name:file:path')
parser.add_argument('-ddb', metavar='name:path', help='add a fully data driven background')
parser.add_argument('-sys', nargs='?', default=[], help='overall sys', metavar='name:dst:up:down', action='append')
parser.add_argument('-histosys', nargs='?', default=[],help='histo sys',metavar='name:dst:up:down', action='append')
parser.add_argument('-sca', nargs='?', default=[], help='norm factor',	metavar='name:dst:value',  action='append')
parser.add_argument('-ssc', nargs='?', default=[], help='shape factor', metavar='name:dst',	   action='append')
parser.add_argument('-fullmc', action='store_true')
parser.add_argument('-combinelumi', action='store_true')
parser.add_argument('-targetlumi', default='1')
parser.add_argument('name')
parser.add_argument('output')
args = parser.parse_args()
sys.argv=[]

name = args.name
targetlumi = args.targetlumi

print "Parse data %s" % args.obs
data={}
data["sname"]= args.obs.split(":")[0]
data["sfile"]= args.obs.split(":")[1]
data["spath"]= args.obs.split(":")[2]
data["hpath"]= dirname(data["spath"])
data["hname"]= os.path.basename(data["spath"])

samples = []
rawoverallsys = []
rawhistosys = []
rawshapesys = []
rawnormfactor = []
rawshapefactor = []
for s in args.sys:
	print "Parse systematic %s" % s
	dic = {}
	dic["name"]=s.split(":")[0]
	dic["dest"]=s.split(":")[1]
	dic["neg"]=float(s.split(":")[2])
	dic["pos"]=float(s.split(":")[3])
	dic["const"]="False"
	rawoverallsys.append(dic)
for s in args.histosys:
	print "Parse histo systematic %s" % s
	dic = {}
	dic["name"]=s.split(":")[0]
	dic["dest"]=s.split(":")[1]
	dic["DN"]=s.split(":")[2]
	dic["UP"]=s.split(":")[3]
	dic["const"]="False"
	rawhistosys.append(dic)
for s in args.sca:
	print "Parse systematic %s" % s
	dic = {}
	dic["name"]=s.split(":")[0]
	dic["dest"]=s.split(":")[1]
	dic["val"]=float(s.split(":")[2])
	dic["neg"]=dic["val"]-0.01
	dic["pos"]=dic["val"]+0.01
	dic["const"]="True"
	rawnormfactor.append(dic)
for s in args.ssc:
	print "Parse systematic %s" % s
	dic = {}
	dic["name"]=s.split(":")[0]
	dic["dest"]=s.split(":")[1]
	dic["const"]="False"
	rawshapefactor.append(dic)
for s in args.sig:
	print "Parse signal %s" % s
	dic = {}
	dic["sname"]=s.split(":")[0]
	dic["spath"]=s.split(":")[1]
	dic["hpath"]=dirname(dic["spath"])
	dic["hname"]=os.path.basename(dic["spath"])
	dic["signal"]=True
	dic["ndd"]="True"
	dic["histosys"] = []
	dic["overallsys"]=[]
	dic["shapesys"]=[]
	dic["normfactor"]=[]
	dic["shapefactor"]=[]
	for sys in rawhistosys:
		if sys["dest"]==dic["sname"]: dic["histosys"].append(sys)
	for sys in rawoverallsys:
		if sys["dest"]==dic["sname"]: dic["overallsys"].append(sys)
	for sys in rawnormfactor:
		if sys["dest"]==dic["sname"]: dic["normfactor"].append(sys)
	for sys in rawshapefactor:
		if sys["dest"]==dic["sname"]: dic["normfactor"].append(sys)
	samples.append(dic)
for s in args.wrsig:
	print "Parse wrong signal %s" % s
	dic = {}
	dic["sname"]=s.split(":")[0]
	dic["spath"]=s.split(":")[1]
	dic["hpath"]=dirname(dic["spath"])
	dic["hname"]=os.path.basename(dic["spath"])
	dic["signal"]=True
	dic["ndd"]="True"
	dic["histosys"] = []
	dic["overallsys"]=[]
	dic["shapesys"]=[]
	dic["normfactor"]=[]
	dic["shapefactor"]=[]
	for sys in rawhistosys:
		if sys["dest"]==dic["sname"]: dic["histosys"].append(sys)
	for sys in rawoverallsys:
		if sys["dest"]==dic["sname"]: dic["overallsys"].append(sys)
	for sys in rawnormfactor:
		if sys["dest"]==dic["sname"]: dic["normfactor"].append(sys)
	for sys in rawshapefactor:
		if sys["dest"]==dic["sname"]: dic["normfactor"].append(sys)
	samples.append(dic)
for s in args.bkg:
	print "Parse background %s" % s
	dic = {}
	dic["sname"]=s.split(":")[0]
	dic["spath"]=s.split(":")[1]
	dic["hpath"]=dirname(dic["spath"])
	dic["hname"]=os.path.basename(dic["spath"])
	dic["signal"]=False
	dic["ndd"]="True"
	dic["histosys"] = []
	dic["overallsys"]=[]
	dic["shapesys"]=[]
	dic["normfactor"]=[]
	dic["shapefactor"]=[]
	for sys in rawoverallsys:
		if sys["dest"]==dic["sname"]: dic["overallsys"].append(sys)
	for sys in rawnormfactor:
		if sys["dest"]==dic["sname"]: dic["normfactor"].append(sys)
	for sys in rawshapefactor:
		if sys["dest"]==dic["sname"]: dic["normfactor"].append(sys)
	for sys in rawhistosys:
		if sys["dest"]==dic["sname"]: dic["histosys"].append(sys)
	samples.append(dic)
if args.ddb:
	print "Create a data driven background"
	fr=ROOT.TFile.Open(data["sfile"],"READ")
	ave_path=data["spath"]
	hname1 = ave_path.replace("data_EM","average")
	hname2 = hname1.replace("data_ME","average")
	hd=fr.Get(hname2)
	bname=args.ddb.split(":")[0]
	bpath=args.ddb.split(":")[1]
	bfile = data["spath"]
	hpath=dirname(bpath)
	hname=os.path.basename(bpath)
	fw=ROOT.TFile(bfile,"UPDATE")
	if hpath:
		fw.mkdir(hpath)
		fw.cd(hpath)
	hb = ROOT.TH1F("%s"%(hname),"%s"%(hname),hd.GetNbinsX(),hd.GetXaxis().GetXmin(),hd.GetXaxis().GetXmax())
	for i in xrange(hd.GetNbinsX()):
		hb.Fill(hd.GetBinCenter(i+1),hd.GetBinContent(i+1))
	hb.Write()
	fw.Close()
	dic = {}
	dic["sname"]="%s"%(hname)
	dic["spath"]="%s/%s"%(hpath,hname) if hpath else "%s"%(hname)
	dic["hpath"]=dirname(dic["spath"])
	dic["hname"]=os.path.basename(dic["spath"])
	dic["signal"]=False
	dic["ndd"]="False"
	dic["shapefactor"]=[{"name":"%s"%(bname),"val":1,"neg":-10,"pos":10,"const":"False"},]
	dic["overallsys"]=[]
	dic["shapesys"]=[]
	dic["histosys"] = []
	dic["normfactor"]=[]
	dic["shapefactor"]=[]
	for sys in rawoverallsys:
		if sys["dest"]==bname: dic["overallsys"].append(sys)
	for sys in rawnormfactor:
		if sys["dest"]==bname: dic["normfactor"].append(sys)
	for sys in rawshapefactor:
		if sys["dest"]==bname: dic["shapefactor"].append(sys)
	for sys in rawhistosys:
		if sys["dest"]==dic["sname"]: dic["histosys"].append(sys)
	samples.append(dic)

#print
if len(os.path.dirname(args.output))>0:
	if not os.path.exists(os.path.dirname(args.output)):
		print "Create output directory: %s" % os.path.dirname(args.output)
		os.mkdir(os.path.dirname(args.output))

fw = open(args.output,"w+")
fw.write("<!DOCTYPE Channel SYSTEM 'HistFactorySchema.dtd'>\n")
fw.write('<Channel Name="%s" InputFile="%s" >\n'%(args.name,data["sfile"]))

fw.write('  <Data InputFile="%(sfile)s" HistoPath="%(hpath)s" HistoName="%(hname)s"/>\n'%data)

if args.fullmc:
    fw.write('	<StatErrorConfig RelErrorThreshold="0.00005" ConstraintType="Poisson"/>\n')

for dic in samples:
	fw.write('  <Sample Name="%(sname)s" NormalizeByTheory="%(ndd)s" HistoPath="%(hpath)s" HistoName="%(hname)s">\n'%dic)
	if args.fullmc:
		fw.write('    <StatError Activate="%(mcstat)s"/>\n'%dic)
	for nui in dic["overallsys"]:
		fw.write('    <OverallSys Name="%(name)s" Low="%(neg).4f" High="%(pos).4f"/>\n' % nui)
	for nui in dic["normfactor"]:
		fw.write('    <NormFactor Name="%(name)s" Val="%(val).4f" Low="%(neg).4f" High="%(pos).4f" Const="%(const)s"/>\n' % nui)
	for nui in dic["histosys"]:
		fw.write('    <HistoSys Name="%(name)s" HistoNameHigh="%(UP)s" HistoNameLow="%(DN)s"/>\n' %nui)
	for nui in dic["shapefactor"]:
		fw.write('    <ShapeFactor Name="%(name)s"/>\n' % nui) 
#Const="%(const)s"/>\n' % nui)
	if dic["signal"]==True:
		fw.write('    <NormFactor Name="mu_BR_htm" Val="0" Low="0" High="200" />\n')
	fw.write('  </Sample>\n')
	
fw.write('</Channel>\n')

print "Generated file: %s" % args.output

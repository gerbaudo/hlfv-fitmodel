#!/usr/bin/env python

import os
import sys

def GetSchema(output):
    output=os.path.dirname(output)
    if not os.path.exists(output):
        os.mkdir(output)
    base=os.path.dirname(output)
    dirs=["/opt/local/etc/root5",]
    if "ROOTSYS" in os.environ:
        dirs.append(os.environ["ROOTSYS"]+"/etc")
    for d in dirs:
        src=d+"/HistFactorySchema.dtd"
        if os.path.exists(src):
            print "Copy %s to %s" % (src,output)
            os.system("cp %s %s" % (src,output))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: %s destination" % os.path.basename(sys.argv[0])
        exit()
    GetSchema(sys.argv[1])

import os
from datetime import datetime
def log(text):
    try:
        pid = os.getpid()
        with open("all_log.txt", 'a') as f:
            f.write(text + " , "+str(pid)+" , " +str(datetime.now())+ " \n")
    except Exception, e:
        print "log write error "+str(e) +"\n"
        
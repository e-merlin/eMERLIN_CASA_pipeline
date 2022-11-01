import time

datapath = ["~/.casa/data"]
rundata = "~/.casa/data/rsync"
logfile = 'casalog-%s.log' % time.strftime("%Y%m%d-%H", time.localtime())
telemetry_enabled = True
crashreporter_enabled = True
nologfile = False
log2term = True
nologger = True
nogui = False
colors = "LightBG"
agg = False
pipeline = False
iplog = True
user_site = False
telemetry_log_directory = "/tmp"
telemetry_log_limit = 20000
telemetry_log_size_interval = 60
telemetry_submit_interval = 604800

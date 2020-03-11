import numpy as np

def rbg_generator(tigId): 
    if tigId is None: return "black"
    np.random.seed((hash(tigId) % 2**30))
    rands = [np.random.rand(), np.random.rand(), np.random.rand()]
    rounded = [round(1e7*r) % 256 for r in rands]
    rgb = ",".join([str(c) for c in rounded])
    return rgb

class FileLogger(object):
    __instance = None
    lines = {}
    flushWhen = 255
    outdir = ""

    def __new__(cls, clean=False, outdir="."):
        if FileLogger.__instance is None or clean:
            FileLogger.__instance = object.__new__(cls)
            FileLogger.__instance.outdir = outdir
            FileLogger.__instance.lines = {}
        return FileLogger.__instance
    
    def flush(self, file):
        with open(self.outdir + "/" + file, "a") as writer:
            for line in self.lines[file]:
                writer.write(line + '\n')
            self.lines[file] = []

    def flush_all(self):
        for file in self.lines:
            if len(self.lines[file]) > 0:
                self.flush(file)

    def write(self, file, line):
        if file not in self.lines:
            self.lines[file] = [line]
        else:
           self.lines[file].append(line)
           if len(self.lines[file]) >= self.flushWhen:
               self.flush(file)
           
    def write_cols(self, file, cols, sep='\t'):
        self.write(file, sep.join([str(c) for c in cols]))
        

def out(text, verboseLevel, param, wait=False):
    if param.VERBOSE >= verboseLevel:
        print(text)
        if wait and param.WAIT:
            input()


NFIX_ATTEMPT = 'nfix'
ALIGNMENT_ATTEMPT = 'alnattempt'
SCAFFOLD_ATTEMPT = 'scaffattempt'

END_ANCHOR_ATTEMPT = 'endanchor'
BLOCK_FILL_ATTEMPT = 'blockfill'
PATH_CLEAN_ATTEMPT = 'pathclean'

SIDE = 'side'
LEFT = 'left'
RIGHT = 'right'

PATH = 'path'
NEXT_PATH = 'nextpath'
SIDE = 'side'
NEXT_SIDE = 'nextside'


REASON = 'reason'
MISSING = 'missing'
ALIGNMENT = 'alignment'
INVALID = 'invalid'
INVALID_INDEX = 'index'
INVALID_ID = 'id'
OVERLAP = 'overlap'
INVALID_STRAND = 'strand'
UNKNOWN = 'unknown'
JOIN_FAIL = 'joinfail'

FLIPPED = 'flipped'

POS = 'pos'
QPOS = 'qpos'
RPOS = 'rpos'

QID = 'qid'
RID = 'rid'

class Report:
    def __init__(self, reportType):
        self.success = True
        self.type = reportType
        self.details = dict()
    
    def add_detail(self, detail, value):
        self.details[detail] = value

    def add_detail_list(self, detail, value):
        self.details.setdefault("detail", []).append(value)

    def set_fail(self, reason=None):
        self.success = False
        if reason is not None:
            self.details[REASON] = reason

    def __repr__(self):
        return self.type + "(" + ("pass" if self.success else "fail") + ")"

    def __str__(self):
        return self.type + "(" + ("pass" if self.success else "fail") + ")\n" + str(self.details)

def print_fail(reports):
    for report in reports:
        if not report.success:
            print(report)
            
class ReportSet:
    def __init__(self, reportType, reports=[]):
        self.reports = reports
        self.type = reportType
        self.details = dict()
    
    def add_detail(self, detail, value):
        self.details[detail] = value

    def success(self):
        for report in self.reports:
            if not report.success:
                return False
        return True
    
    def has_success(self):
        for report in self.reports:
            if report.success:
                return True
        return False

    
    def __getitem__(self, key):
         return self.reports[key]
    def __repr__(self):
        return self.type + " (" + str(sum([1 if report.success else 0 for report in self.reports])) + "/" + str(len(self.reports)) + ")"
    def __str__(self):
        return self.type + " (" + str(sum([1 if report.success else 0 for report in self.reports])) + "/" + str(len(self.reports)) + ")"

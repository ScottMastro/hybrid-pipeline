class Path:
    def __init__(self):
        self.path = []
        
    def add_fork(self, fork):
        self.path.append(fork)
    
    def add_path(self, other):
        self.path.extend(other.path)
    
    def last_fork(self):
        if len(self.path) == 0:
            return None
        return self.path[-1]

    def __repr__(self):
        return "\n".join([f.__repr__() for f in self.path])
    def __str__(self):
        return "\n".join([f.__str__() for f in self.path])


class Fork:
    def __init__(self, qid, qpos, qstrand, rid, rpos, rstrand):
        self.qid = qid     
        self.qpos = qpos     
        self.qstrand = qstrand
        
        self.rid = rid     
        self.rpos = rpos
        self.rstrand = rstrand
        self.switch = "q"
                
    def switch_query(self):
        self.switch = "q"

    def switch_reference(self):
        self.switch = "r"

    def get_pos(self, q=True):
        if q: return self.qpos
        else: return self.rpos
        
    def __repr__(self):
        return str(self.qid) + " " + str(self.qpos) + " (" + str(self.qstrand) + ") " + \
            str(self.rid) + " " + str(self.rpos) + " (" + str(self.rstrand) + ")" 
    def __str__(self):
        return self.__repr__()
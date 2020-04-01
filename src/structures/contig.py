class Contig:

    def __init__(self, id, size, mblocks=[]):
        self.id=id
        self.size=size
        self.mblocks = mblocks

    def is_empty(self): return len(self.mblocks) == 0
   
    def __len__(self): return len(self.mblocks)
    def __repr__(self):
        return str(self.mblocks) 
    def __str__(self):
        return "contig id=" + str(self.id)


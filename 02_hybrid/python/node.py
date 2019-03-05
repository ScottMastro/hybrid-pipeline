class Node:
    id=""
    left_used=False
    right_used=False
    left=dict()
    right=dict()

    def __init__(self, tigId):
        self.id=tigId

    def get_right(self):
        if not self.right_used:
            self.right_used=True
            return self.right
        else:
            return dict()
    def get_left(self):
        if not self.left_used:
            self.left_used=True
            return self.left
        else:
            return dict()

    def get_id(self):
        return self.id
    def set_used(self):
        self.used=True
    def is_used(self):
        return self.left_used or self.right_used
    
    def get_next(self, prevId):
        
        if prevId in self.right:
            self.right_used = True
            return self.get_left()
        if prevId in self.left:
            self.left_used = True
            return self.get_right()
        
        if prevId == "__left__":
            return self.get_left()
        if prevId == "__right__":
            return self.get_right()
        
        print("WARNING: prevId not found")

    
    #hint=0, left contig / hint=1, right contig
    def add_edge(self, bridge, hint=0):
        b0, b1 = bridge
        
        if hint == 0:
            dir = b0[3]
            nextId = b1[0]
            if dir == 1:
                self.right[nextId] = bridge
            if dir == -1:
                self.left[nextId] = bridge
            
        if hint == 1:
            dir = b1[3]
            nextId = b0[0]
            if dir == 1:
                self.left[nextId] = bridge
            if dir == -1:
                self.right[nextId] = bridge
            
    def is_singleton(self):
        return (len(self.left) + len(self.right)) == 0
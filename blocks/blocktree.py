class Node:
    data=None
    def __init__(self, data):

        self.left = None
        self.right = None
        self.data = data

    def insert(self, node):
        if node < self:
            if self.left is None:
                self.left = node
            else:
                self.left.insert(node)
        elif not node < self:
            if self.right is None:
                self.right = node
            else:
                self.right.insert(node)

    def __lt__(self, other):
        return self.data[0] < other.data[0]

    def __str__(self):
        string = "left:\n"
        if self.left is not None:
            string = string + self.left.__repr__()
        string = string + "\nroot:\n" + self.data.__repr__() + "\nright:\n"
        if self.right is not None:
            string = string + self.right.__repr__()
        return string
    
    def __repr__(self):
        string = ""
        if self.left is not None:
            string = string + self.left.__repr__()
        string = string + self.data.__repr__()
        if self.right is not None:
            string = string + self.right.__repr__()
        return string

def traverse(node, prevNodeL, prevNodeR, collapseFn):

    keep, trash = [], []
    flagL, flagR = True, True
    
    if prevNodeL is not None:
        flagL = collapseFn(prevNodeL, node)
    if prevNodeR is not None:
        flagR = collapseFn(node, prevNodeR)

    if flagL and flagR:
        nextNodeL = node
        nextNodeR = node
        keep.append(node.data[-1])
    else:
        nextNodeL = prevNodeL
        nextNodeR = prevNodeR
        trash.append(node.data[-1])
            
    left = right = []
    leftTrash = rightTrash = []

    if node.left is not None:
        left, leftTrash = traverse(node.left, prevNodeL, nextNodeR, collapseFn)
    if node.right is not None:
        right, rightTrash = traverse(node.right, nextNodeL, prevNodeR, collapseFn)
        
    return (left + keep + right, leftTrash + trash + rightTrash)
    
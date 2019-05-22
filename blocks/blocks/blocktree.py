def display(node):
    lines, _, _, _ = display_aux(node)
    string = "\n".join([line for line in lines])
    return string

def display_aux(node):
    """Returns list of strings, width, height, and horizontal coordinate of the root."""
    # No child.
    if node.r is None and node.l is None:
        line = node.data.__repr__()
        width = len(line)
        height = 1
        middle = width // 2
        return [line], width, height, middle

    # Only left child.
    if node.r is None:
        lines, n, p, x = display_aux(node.l)
        s = node.data.__repr__()
        u = len(s)
        first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s
        second_line = x * ' ' + '/' + (n - x - 1 + u) * ' '
        shifted_lines = [line + u * ' ' for line in lines]
        return [first_line, second_line] + shifted_lines, n + u, p + 2, n + u // 2

    # Only right child.
    if node.l is None:
        lines, n, p, x = display_aux(node.r)
        s = node.data.__repr__()
        u = len(s)
        first_line = s + x * '_' + (n - x) * ' '
        second_line = (u + x) * ' ' + '\\' + (n - x - 1) * ' '
        shifted_lines = [u * ' ' + line for line in lines]
        return [first_line, second_line] + shifted_lines, n + u, p + 2, u // 2

    # Two children.
    left, n, p, x = display_aux(node.l)
    right, m, q, y = display_aux(node.r)
    s = node.data.__repr__()
    u = len(s)
    first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s + y * '_' + (m - y) * ' '
    second_line = x * ' ' + '/' + (n - x - 1 + u + y) * ' ' + '\\' + (m - y - 1) * ' '
    if p < q:
        left += [n * ' '] * (q - p)
    elif q < p:
        right += [m * ' '] * (p - q)
    zipped_lines = zip(left, right)
    lines = [first_line, second_line] + [a + u * ' ' + b for a, b in zipped_lines]
    return lines, n + m + u, max(p, q) + 2, n + u // 2



class Node:
    def __init__(self, data, q=True):
        self.data = data
        self.l = None
        self.r = None
        self.q = q
        
    def insert(self, node):
        if node < self:
            if self.l is None:
                self.l = node
            else:
                self.l.insert(node)
        elif not node < self:
            if self.r is None:
                self.r = node
            else:
                self.r.insert(node)

    def left_pos(self, q=True):
        return self.data.left(q)
    def right_pos(self, q=True): 
        return self.data.right(q)
    def size(self, q=True): 
        return self.data.sum(self.q)
    def span(self, q=True): 
        return self.data.span(self.q)

    def to_list(self):
        return self.l.to_list() if self.l is not None else [] + \
                [self.data] + \
                self.r.to_list() if self.r is not None else []
        
    def __lt__(self, other):
        return self.left_pos(self.q) < other.left_pos(self.q)

    def __str__(self):
        return display(self)
    
    def __repr__(self):
        return self.data.__repr__()
    
'''
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
'''

def traverse(node, prevNodeL, prevNodeR, collapseFn):

    if node is None:
        return([],[])
    
    flagL, flagR = True, True
    flagLrecurse, flagRrecurse = False, False

    if prevNodeL is not None:
        flagL = collapseFn(prevNodeL, node)
        
        #continue down tree to try and recover node
        if not flagL:  
            left, leftTrash = traverse(node.l, prevNodeL, node, collapseFn)
            flagLrecurse = True
            if len(left) > 0:
                flagL = collapseFn(left[-1], node)    
        
    if prevNodeR is not None:
        flagR = collapseFn(node, prevNodeR)    
    
        #continue down tree to try and recover node
        if not flagR:
            right, rightTrash = traverse(node.r, node, prevNodeR, collapseFn)
            flagRrecurse = True
            if len(right) > 0:
                flagR = collapseFn(node, right[0])    

    if flagL and flagR:
        if not flagLrecurse:
            left, leftTrash = traverse(node.l, prevNodeL, node, collapseFn)
        if not flagRrecurse:
            right, rightTrash = traverse(node.r, node, prevNodeR, collapseFn)
        return (left + [node] + right, leftTrash + rightTrash)
    else:
        left, leftTrash = traverse(node.l, prevNodeL, prevNodeR, collapseFn)
        right, rightTrash = traverse(node.r, prevNodeL, prevNodeR, collapseFn)
        return (left + right, leftTrash + [node] + rightTrash)

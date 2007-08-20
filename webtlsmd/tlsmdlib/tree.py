## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

class Tree(object):
    def __init__(self, name = None):
        self.__name = name
        self.__parent = None 
        self.__children = []

    def __str__(self):
        return str(self.__name)

    def empty(self):
        return len(self.__children) == 0

    def depth(self):
        depth = 0
        for d, c in self.iter_depth_first():
            depth = max(depth, d)
        return depth

    def append(self, child):
        child.__parent = self
        self.__children.append(child)

    def unparent(self):
        self.__parent = None

    def parent(self):
        return self.__parent

    def iter_children(self):
        return iter(self.__children)

    def iter_depth_first(self):
        depth = 1
        children = self.__children
        while len(children) > 0:
            q = list()
            for child in children:
                yield depth, child
                for qchild in child.__children:
                    q.append(qchild)
            children = q
            depth += 1

    def iter_depth(self, depth):
        ic = self.iter_children()
        d, c = ic.next()
        while d < depth:
            d, c = ic.next()
        while d == depth:
            yield c
            d, c = ic.next()

def testmain():
    def add_tree(tree, num_children, depth):
        if depth < 1:
            return
        for i in xrange(num_children):
            child = Tree(depth)
            tree.append(child)
            add_tree(child, num_children, depth - 1)

    root = Tree()
    add_tree(root, 2, 4)

    for depth, tree in root.iter_children():
        print "%s:%d " % (tree, depth),
    print
    for tree in root.iter_depth(3):
        print tree,
    print

if __name__ == "__main__":
    testmain()

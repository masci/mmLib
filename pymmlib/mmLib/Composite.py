from __future__ import generators

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.


CompositeError = 'CompositeError'



class Composite:
    """Base class for all containers used in mmLib.Structure.  This implements
    a Composite design pattern tree of containers."""

    def __init__(self):
        self.__parent     = None
        self.__collapsed  = 0
        self.__child_list = []


    def __str__(self):
        return "Composite"


    def blockingIterator(self, block_func = None):
        """Preorder Traversal for a Composite tree using Python generators."""
        yield self

        if self.__collapsed:
            return

        if block_func and block_func(self):
            return
        
        for child in self.__child_list:
            for x in child.blockingIterator(block_func):
                yield x


    def __iter__(self):
        """Preorder Traversal for a Composite tree using Python generators."""
        yield self

        if self.__collapsed:
            return
        
        for child in self.__child_list:
            for x in child: yield x


    def isCollapsed(self):
        return self.__collapsed


    def Collapse(self):
        self.__collapsed = 1


    def Expand(self):
        self.__collapsed = 0


    def getDepth(self):
        """Returns the depth, the root composite is depth 0."""
        depth = 0
        ancestor = self.__parent
        while ancestor:
            depth += 1
            ancestor = ancestor.parent
        return depth


    def getDegree(self):
        """Returns the number of children (degree)."""
        return len(self.__child_list)


    def countDescendants(self):
        """Counts all decendants."""
        n = self.getDegree()
        for child in self.__child_list:
            n += child.countDescendants()
        return n


    def getRoot(self):
        """Returns the root composite."""
        composite = self
        while composite.__parent:
            composite = composite.__parent
        return composite


    def getParent(self):
        return self.__parent


    def getChildList(self):
        return self.__child_list


    def getPath(self):
        """Returns the tree-path to the composite as a list of its
        parent composites."""
        path_list = [self]
        parent = self.__parent
        while parent:
            path_list.insert(0, parent)
            parent = parent.__parent
        return path_list


    def getIndexPath(self):
        """Returns the tree-path to the composite as a list of its
        integer indexes."""
        ipath_list = []
        child = self
        parent = child.__parent
        while parent:
            ipath_list.insert(0, parent.__child_list.index(child))
            child = parent
            parent = parent.__parent
        return ipath_list


    def getParentList(self):
        """Returns a list of the parent composites back to the root."""
        list = []
        composite = self
        while composite.__parent:
            composite = composite.__parent
            list.append(composite)
        return list


    def getLowestCommonAncestor(self, composite):
        """Returns the lowest common ancesotry of self and argument
        composite."""
        pl1 = self.getParentList()
        pl2 = composite.getParentList()
        pl1.reverse()
        pl2.reverse()

        ancestor = None
        for i in range(min(len(pl1), len(pl2))):
            if pl1[i] == pl2[i]:
                ancestor = pl1[i]
            else:
                break

        return ancestor

        
    def isDescendantOf(self, composite):
        """Returns true if self composite is a decent of argument composite."""
        ancestor = self.__parent
        while ancestor:
            if ancestor == composite:
                return 1
            ancestor = ancestor.__parent
        return 0


    def removeChild(self, child):
        """Removes child, ignores errors."""
        try:
            self.__child_list.remove(child)
        except ValueError:
            return
        child.parent = None


    def prependChild(self, child):
        if child == self:
            raise CompositeError, 'Prepending self not allowed'
        if self.isDescendantOf(child):
            raise CompositeError, 'Prepending of parents not allowed'

        if child.parent:
            child.__parent.removeChild(child)
        child.__parent = self

        self.__child_list.insert(0, child)

            
    def appendChild(self, child):
        if child == self:
            raise CompositeError, 'Appending self not allowed'
        if self.isDescendantOf(child):
            raise CompositeError, 'Appending of parents not allowed'

        if child.__parent:
            child.__parent.removeChild(child)
        child.__parent = self
            
        self.__child_list.append(child)


    def insertBefore(self, composite):
        """Inserts composite as the sibling before self composite."""
        if composite == self:
            raise CompositeError, 'Insert self not allowed'
        if self.isDescendantOf(composite):
            raise CompositeError, 'Insert of parents not allowed'
        if not self.__parent:
            raise CompositeError, 'Cannot insert sibling on root composite'

        if composite.__parent:
            composite.__parent.removeChild(composite)
        composite.__parent = self.__parent

        i = self.__parent.__child_list.index(self)
        self.__parent.__child_list.insert(i - 1, composite)


    def insertAfter(self, composite):
        """Inserts composite as the sibling before self composite."""
        if composite == self:
            raise CompositeError, 'Insert self not allowed'
        if self.isDescendantOf(composite):
            raise CompositeError, 'Insert of parents not allowed'
        if not self.__parent:
            raise CompositeError, 'Cannot insert sibling on root composite'

        if composite.__parent:
            composite.__parent.removeChild(composite)
        composite.__parent = self.__parent

        i = self.__parent.__child_list.index(self)
        self.__parent.__child_list.insert(i, composite)


    def getParent(self):
        return self.__parent


    def insertParent(self, parent, first, last):
        if not first.__parent:
            raise CompositeError, 'First composite has no parent'
        if not last.__parent:
            raise CompositeError, 'Last composite has no parent'
        if first.parent != last.parent:
            raise CompositeError, 'Parents of first and last do not match'
        if first.isDescendantOf(parent):
            raise CompositeError, 'First/Last already descents of parent'

        self.__parent = first.__parent
        fi = self.__parent.__child_list.index(first)
        li = self.__parent.__child_list.index(last)

        self.__child_list = self.__parent.__child_list[fi:fi+1]
        self.__parent.__child_list[fi:fi+1] = [self]

        for child in self.__child_list:
            child.__parent = self


    def getFirstChild(self):
        try:
            return self.__child_list[0]
        except IndexError:
            raise IndexError, 'Composite has no children'


    def getLastChild(self):
        try:
            return self.__child_list[-1]
        except IndexError:
            raise IndexError, 'Composite has no children'


    def getChild(self, index):
        try:
            return self.__child_list[index]
        except IndexError:
            raise IndexError, 'No child at given index=%d' % (index)


    def getIndex(self, child):
        return self.__child_list.index(child)


    def getSibling(self, delta_index):
        if not self.__parent:
            raise IndexError, 'Root composite has no siblings'
        index = self.__parent.__child_list.index(self) + delta_index
        return self.__parent.getChild(index)



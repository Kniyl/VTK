"""Test vtkVariant support in VTK-Python

The following vtkVariant features have to be tested:
- Construction from various types
- Automatic arg conversion for methods with vtkVariant args
- Access of various types
- Operators < <= == > >=
- Use of vtkVariant as a dict key

The following features are not supported
- The ToNumeric method

Created on May 12, 2010 by David Gobbi

"""

import sys
import exceptions
import vtk
from vtk.test import Testing

unicode_support = False
try:
    unicode('hello')
    unicode_support = True
except:
    print "unicode not supported on this python installation"


class TestVariant(Testing.vtkTest):
    def testDefaultConstructor(self):
        """Default constructor"""
        v = vtk.vtkVariant()
        self.assertEqual(v.IsValid(), False)
        self.assertEqual(v.GetType(), 0)

    def testCopyConstructor(self):
        """Construct from another vtkVariant"""
        u = vtk.vtkVariant('test')
        v = vtk.vtkVariant(u)
        self.assertEqual(v.GetType(), vtk.VTK_STRING)
        self.assertEqual(v.ToString(), u.ToString())

    def testIntConstructor(self):
        """Construct from int"""
        v = vtk.vtkVariant(10)
        self.assertEqual(v.GetType(), vtk.VTK_INT)
        self.assertEqual(v.ToInt(), 10)

    def testFloatConstructor(self):
        """Construct from float"""
        v = vtk.vtkVariant(10.0)
        self.assertEqual(v.GetType(), vtk.VTK_DOUBLE)
        self.assertEqual(v.ToDouble(), 10.0)

    def testStringConstructor(self):
        """Construct from string"""
        v = vtk.vtkVariant('hello')
        self.assertEqual(v.GetType(), vtk.VTK_STRING)
        self.assertEqual(v.ToString(), 'hello')

    def testUnicodeConstructor(self):
        """Construct from unicode"""
        if unicode_support:
            v = vtk.vtkVariant(u'hello')
            self.assertEqual(v.GetType(), vtk.VTK_UNICODE_STRING)
            self.assertEqual(v.ToUnicodeString(), u'hello')

    def testObjectConstructor(self):
        """Construct from VTK object"""
        o = vtk.vtkIntArray()
        v = vtk.vtkVariant(o)
        self.assertEqual(v.GetType(), vtk.VTK_OBJECT)
        self.assertEqual(v.GetTypeAsString(), o.GetClassName())
        self.assertEqual(v.ToVTKObject(), o)

    def testTwoArgConstructor(self):
        """Construct with a specific type"""
        # construct with conversion to int
        v = vtk.vtkVariant('10')
        u = vtk.vtkVariant(v, vtk.VTK_INT)
        self.assertEqual(u.GetType(), vtk.VTK_INT)
        self.assertEqual(v.ToInt(), u.ToInt())

        # construct with conversion to double
        v = vtk.vtkVariant(10)
        u = vtk.vtkVariant(v, vtk.VTK_DOUBLE)
        self.assertEqual(u.GetType(), vtk.VTK_DOUBLE)
        self.assertEqual(u.ToDouble(), 10.0)

        # failed conversion to vtkObject
        v = vtk.vtkVariant(10)
        u = vtk.vtkVariant(v, vtk.VTK_OBJECT)
        self.assertEqual(u.IsValid(), False)

    def testAutomaticArgConversion(self):
        """Automatic construction of variants to resolve args"""
        # use with one of vtkVariant's own constructors
        v = vtk.vtkVariant('10', vtk.VTK_INT)
        self.assertEqual(v.ToInt(), 10)
        self.assertEqual(v.GetType(), vtk.VTK_INT)

        # use with vtkVariantArray
        a = vtk.vtkVariantArray()
        i = a.InsertNextValue(10)
        v = a.GetValue(i)
        self.assertEqual(v.GetType(), vtk.VTK_INT)
        self.assertEqual(v.ToInt(), 10)
        i = a.InsertNextValue(10.0)
        v = a.GetValue(i)
        self.assertEqual(v.GetType(), vtk.VTK_DOUBLE)
        self.assertEqual(v.ToDouble(), 10.0)
        i = a.InsertNextValue('10')
        v = a.GetValue(i)
        self.assertEqual(v.GetType(), vtk.VTK_STRING)
        self.assertEqual(v.ToString(), '10')

    def testCompare(self):
        """Use comparison operators to sort an array of vtkVariants"""
        if not unicode_support:
            return
        l = map(vtk.vtkVariant, [1, 2.5, vtk.vtkVariant(), "0", u'hello'])
        s = map(vtk.vtkVariant, [vtk.vtkVariant(), "0", 1, 2.5, u'hello'])
        l.sort()
        self.assertEqual(l, s)

    def testVariantExtract(self):
        """Use vtkVariantExtract"""
        l = [1, '2', u'3', 4.0, vtk.vtkObject()]
        s = map(vtk.vtkVariant, l)
        m = map(vtk.vtkVariantExtract, s)
        self.assertEqual(l, m)

    def testHash(self):
        """Use a variant as a dict key"""
        d = {}
        # doubles, ints, srings, all hash as strings
        d[vtk.vtkVariant(1.0)] = 'double'
        d[vtk.vtkVariant(1)] = 'int'
        self.assertEqual(d[vtk.vtkVariant('1')], 'int')

        # strings and unicode have the same hash
        d[vtk.vtkVariant(u's')] = 'unicode'
        d[vtk.vtkVariant('s')] = 'string'
        self.assertEqual(d[vtk.vtkVariant(u's')], 'string')

        # every vtkObject is hashed by memory address
        o1 = vtk.vtkIntArray()
        o2 = vtk.vtkIntArray()
        d[vtk.vtkVariant(o1)] = 'vtkIntArray1'
        d[vtk.vtkVariant(o2)] = 'vtkIntArray2'
        self.assertEqual(d[vtk.vtkVariant(o1)], 'vtkIntArray1')
        self.assertEqual(d[vtk.vtkVariant(o2)], 'vtkIntArray2')

        # invalid variants all hash the same
        d[vtk.vtkVariant()] = 'invalid'
        self.assertEqual(d[vtk.vtkVariant()], 'invalid')

    def testPassByValueReturnByReference(self):
        """Pass vtkVariant by value, return by reference"""
        a = vtk.vtkVariantArray()
        a.SetNumberOfValues(1)
        v = vtk.vtkVariant(1)
        a.SetValue(0, v)
        u = a.GetValue(0)
        self.assertEqual(u.ToInt(), v.ToInt())
        self.assertEqual(u.GetType(), v.GetType())
        self.assertEqual(u.IsValid(), v.IsValid())

    def testPassByReferenceReturnByValue(self):
        """Pass vtkVariant by reference, return by value."""
        a = vtk.vtkArray.CreateArray(1, vtk.VTK_INT)
        a.Resize(1,1)
        v = vtk.vtkVariant(1)
        a.SetVariantValue(0, 0, v)
        u = a.GetVariantValue(0, 0)
        self.assertEqual(u.ToInt(), v.ToInt())
        self.assertEqual(u.GetType(), v.GetType())
        self.assertEqual(u.IsValid(), v.IsValid())

if __name__ == "__main__":
    Testing.main([(TestVariant, 'test')])

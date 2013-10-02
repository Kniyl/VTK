# Copyright (c) Twisted Matrix Laboratories.
# See LICENSE for details.

"""
Test cases for L{jelly} object serialization.
"""

import datetime

try:
    import decimal
except ImportError:
    decimal = None

from twisted.spread import jelly, pb
from twisted.python.compat import set, frozenset
from twisted.trial import unittest
from twisted.test.proto_helpers import StringTransport


class TestNode(object, jelly.Jellyable):
    """
    An object to test jellyfying of new style class instances.
    """
    classAttr = 4

    def __init__(self, parent=None):
        if parent:
            self.id = parent.id + 1
            parent.children.append(self)
        else:
            self.id = 1
        self.parent = parent
        self.children = []



class A:
    """
    Dummy class.
    """

    def amethod(self):
        """
        Method tp be used in serialization tests.
        """



def afunc(self):
    """
    A dummy function to test function serialization.
    """



class B:
    """
    Dummy class.
    """

    def bmethod(self):
        """
        Method to be used in serialization tests.
        """



class C:
    """
    Dummy class.
    """

    def cmethod(self):
        """
        Method to be used in serialization tests.
        """



class D(object):
    """
    Dummy new-style class.
    """



class E(object):
    """
    Dummy new-style class with slots.
    """

    __slots__ = ("x", "y")

    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y


    def __getstate__(self):
        return {"x" : self.x, "y" : self.y}


    def __setstate__(self, state):
        self.x = state["x"]
        self.y = state["y"]



class SimpleJellyTest:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def isTheSameAs(self, other):
        return self.__dict__ == other.__dict__



class JellyTestCase(unittest.TestCase):
    """
    Testcases for L{jelly} module serialization.

    @cvar decimalData: serialized version of decimal data, to be used in tests.
    @type decimalData: C{list}
    """

    def _testSecurity(self, inputList, atom):
        """
        Helper test method to test security options for a type.

        @param inputList: a sample input for the type.
        @param inputList: C{list}

        @param atom: atom identifier for the type.
        @type atom: C{str}
        """
        c = jelly.jelly(inputList)
        taster = jelly.SecurityOptions()
        taster.allowBasicTypes()
        # By default, it should succeed
        jelly.unjelly(c, taster)
        taster.allowedTypes.pop(atom)
        # But it should raise an exception when disallowed
        self.assertRaises(jelly.InsecureJelly, jelly.unjelly, c, taster)


    def test_methodSelfIdentity(self):
        a = A()
        b = B()
        a.bmethod = b.bmethod
        b.a = a
        im_ = jelly.unjelly(jelly.jelly(b)).a.bmethod
        self.assertEqual(im_.im_class, im_.im_self.__class__)


    def test_methodsNotSelfIdentity(self):
        """
        If a class change after an instance has been created, L{jelly.unjelly}
        shoud raise a C{TypeError} when trying to unjelly the instance.
        """
        a = A()
        b = B()
        c = C()
        a.bmethod = c.cmethod
        b.a = a
        savecmethod = C.cmethod
        del C.cmethod
        try:
            self.assertRaises(TypeError, jelly.unjelly, jelly.jelly(b))
        finally:
            C.cmethod = savecmethod


    def test_newStyle(self):
        n = D()
        n.x = 1
        n2 = D()
        n.n2 = n2
        n.n3 = n2
        c = jelly.jelly(n)
        m = jelly.unjelly(c)
        self.assertIsInstance(m, D)
        self.assertIdentical(m.n2, m.n3)


    def test_newStyleWithSlots(self):
        """
        A class defined with I{slots} can be jellied and unjellied with the
        values for its attributes preserved.
        """
        n = E()
        n.x = 1
        c = jelly.jelly(n)
        m = jelly.unjelly(c)
        self.assertIsInstance(m, E)
        self.assertEqual(n.x, 1)


    def test_typeOldStyle(self):
        """
        Test that an old style class type can be jellied and unjellied
        to the original type.
        """
        t = [C]
        r = jelly.unjelly(jelly.jelly(t))
        self.assertEqual(t, r)


    def test_typeNewStyle(self):
        """
        Test that a new style class type can be jellied and unjellied
        to the original type.
        """
        t = [D]
        r = jelly.unjelly(jelly.jelly(t))
        self.assertEqual(t, r)


    def test_typeBuiltin(self):
        """
        Test that a builtin type can be jellied and unjellied to the original
        type.
        """
        t = [str]
        r = jelly.unjelly(jelly.jelly(t))
        self.assertEqual(t, r)


    def test_dateTime(self):
        dtn = datetime.datetime.now()
        dtd = datetime.datetime.now() - dtn
        input = [dtn, dtd]
        c = jelly.jelly(input)
        output = jelly.unjelly(c)
        self.assertEqual(input, output)
        self.assertNotIdentical(input, output)


    def test_decimal(self):
        """
        Jellying L{decimal.Decimal} instances and then unjellying the result
        should produce objects which represent the values of the original
        inputs.
        """
        inputList = [decimal.Decimal('9.95'),
                     decimal.Decimal(0),
                     decimal.Decimal(123456),
                     decimal.Decimal('-78.901')]
        c = jelly.jelly(inputList)
        output = jelly.unjelly(c)
        self.assertEqual(inputList, output)
        self.assertNotIdentical(inputList, output)


    decimalData = ['list', ['decimal', 995, -2], ['decimal', 0, 0],
                           ['decimal', 123456, 0], ['decimal', -78901, -3]]


    def test_decimalUnjelly(self):
        """
        Unjellying the s-expressions produced by jelly for L{decimal.Decimal}
        instances should result in L{decimal.Decimal} instances with the values
        represented by the s-expressions.

        This test also verifies that C{self.decimalData} contains valid jellied
        data.  This is important since L{test_decimalMissing} re-uses
        C{self.decimalData} and is expected to be unable to produce
        L{decimal.Decimal} instances even though the s-expression correctly
        represents a list of them.
        """
        expected = [decimal.Decimal('9.95'),
                    decimal.Decimal(0),
                    decimal.Decimal(123456),
                    decimal.Decimal('-78.901')]
        output = jelly.unjelly(self.decimalData)
        self.assertEqual(output, expected)


    def test_decimalMissing(self):
        """
        If decimal is unavailable on the unjelly side, L{jelly.unjelly} should
        gracefully return L{jelly.Unpersistable} objects.
        """
        self.patch(jelly, 'decimal', None)
        output = jelly.unjelly(self.decimalData)
        self.assertEqual(len(output), 4)
        for i in range(4):
            self.assertIsInstance(output[i], jelly.Unpersistable)
        self.assertEqual(output[0].reason,
            "Could not unpersist decimal: 9.95")
        self.assertEqual(output[1].reason,
            "Could not unpersist decimal: 0")
        self.assertEqual(output[2].reason,
            "Could not unpersist decimal: 123456")
        self.assertEqual(output[3].reason,
            "Could not unpersist decimal: -78.901")


    def test_decimalSecurity(self):
        """
        By default, C{decimal} objects should be allowed by
        L{jelly.SecurityOptions}. If not allowed, L{jelly.unjelly} should raise
        L{jelly.InsecureJelly} when trying to unjelly it.
        """
        inputList = [decimal.Decimal('9.95')]
        self._testSecurity(inputList, "decimal")

    if decimal is None:
        skipReason = "decimal not available"
        test_decimal.skip = skipReason
        test_decimalUnjelly.skip = skipReason
        test_decimalSecurity.skip = skipReason


    def test_set(self):
        """
        Jellying C{set} instances and then unjellying the result
        should produce objects which represent the values of the original
        inputs.
        """
        inputList = [set([1, 2, 3])]
        output = jelly.unjelly(jelly.jelly(inputList))
        self.assertEqual(inputList, output)
        self.assertNotIdentical(inputList, output)


    def test_frozenset(self):
        """
        Jellying C{frozenset} instances and then unjellying the result
        should produce objects which represent the values of the original
        inputs.
        """
        inputList = [frozenset([1, 2, 3])]
        output = jelly.unjelly(jelly.jelly(inputList))
        self.assertEqual(inputList, output)
        self.assertNotIdentical(inputList, output)


    def test_setSecurity(self):
        """
        By default, C{set} objects should be allowed by
        L{jelly.SecurityOptions}. If not allowed, L{jelly.unjelly} should raise
        L{jelly.InsecureJelly} when trying to unjelly it.
        """
        inputList = [set([1, 2, 3])]
        self._testSecurity(inputList, "set")


    def test_frozensetSecurity(self):
        """
        By default, C{frozenset} objects should be allowed by
        L{jelly.SecurityOptions}. If not allowed, L{jelly.unjelly} should raise
        L{jelly.InsecureJelly} when trying to unjelly it.
        """
        inputList = [frozenset([1, 2, 3])]
        self._testSecurity(inputList, "frozenset")


    def test_oldSets(self):
        """
        Test jellying C{sets.Set}: it should serialize to the same thing as
        C{set} jelly, and be unjellied as C{set} if available.
        """
        inputList = [jelly._sets.Set([1, 2, 3])]
        inputJelly = jelly.jelly(inputList)
        self.assertEqual(inputJelly, jelly.jelly([set([1, 2, 3])]))
        output = jelly.unjelly(inputJelly)
        # Even if the class is different, it should coerce to the same list
        self.assertEqual(list(inputList[0]), list(output[0]))
        if set is jelly._sets.Set:
            self.assertIsInstance(output[0], jelly._sets.Set)
        else:
            self.assertIsInstance(output[0], set)


    def test_oldImmutableSets(self):
        """
        Test jellying C{sets.ImmutableSet}: it should serialize to the same
        thing as C{frozenset} jelly, and be unjellied as C{frozenset} if
        available.
        """
        inputList = [jelly._sets.ImmutableSet([1, 2, 3])]
        inputJelly = jelly.jelly(inputList)
        self.assertEqual(inputJelly, jelly.jelly([frozenset([1, 2, 3])]))
        output = jelly.unjelly(inputJelly)
        # Even if the class is different, it should coerce to the same list
        self.assertEqual(list(inputList[0]), list(output[0]))
        if frozenset is jelly._sets.ImmutableSet:
            self.assertIsInstance(output[0], jelly._sets.ImmutableSet)
        else:
            self.assertIsInstance(output[0], frozenset)


    def test_simple(self):
        """
        Simplest test case.
        """
        self.failUnless(SimpleJellyTest('a', 'b').isTheSameAs(
                        SimpleJellyTest('a', 'b')))
        a = SimpleJellyTest(1, 2)
        cereal = jelly.jelly(a)
        b = jelly.unjelly(cereal)
        self.failUnless(a.isTheSameAs(b))


    def test_identity(self):
        """
        Test to make sure that objects retain identity properly.
        """
        x = []
        y = (x)
        x.append(y)
        x.append(y)
        self.assertIdentical(x[0], x[1])
        self.assertIdentical(x[0][0], x)
        s = jelly.jelly(x)
        z = jelly.unjelly(s)
        self.assertIdentical(z[0], z[1])
        self.assertIdentical(z[0][0], z)


    def test_unicode(self):
        x = unicode('blah')
        y = jelly.unjelly(jelly.jelly(x))
        self.assertEqual(x, y)
        self.assertEqual(type(x), type(y))


    def test_stressReferences(self):
        reref = []
        toplevelTuple = ({'list': reref}, reref)
        reref.append(toplevelTuple)
        s = jelly.jelly(toplevelTuple)
        z = jelly.unjelly(s)
        self.assertIdentical(z[0]['list'], z[1])
        self.assertIdentical(z[0]['list'][0], z)


    def test_moreReferences(self):
        a = []
        t = (a,)
        a.append((t,))
        s = jelly.jelly(t)
        z = jelly.unjelly(s)
        self.assertIdentical(z[0][0][0], z)


    def test_typeSecurity(self):
        """
        Test for type-level security of serialization.
        """
        taster = jelly.SecurityOptions()
        dct = jelly.jelly({})
        self.assertRaises(jelly.InsecureJelly, jelly.unjelly, dct, taster)


    def test_newStyleClasses(self):
        j = jelly.jelly(D)
        uj = jelly.unjelly(D)
        self.assertIdentical(D, uj)


    def test_lotsaTypes(self):
        """
        Test for all types currently supported in jelly
        """
        a = A()
        jelly.unjelly(jelly.jelly(a))
        jelly.unjelly(jelly.jelly(a.amethod))
        items = [afunc, [1, 2, 3], not bool(1), bool(1), 'test', 20.3,
                 (1, 2, 3), None, A, unittest, {'a': 1}, A.amethod]
        for i in items:
            self.assertEqual(i, jelly.unjelly(jelly.jelly(i)))


    def test_setState(self):
        global TupleState
        class TupleState:
            def __init__(self, other):
                self.other = other
            def __getstate__(self):
                return (self.other,)
            def __setstate__(self, state):
                self.other = state[0]
            def __hash__(self):
                return hash(self.other)
        a = A()
        t1 = TupleState(a)
        t2 = TupleState(a)
        t3 = TupleState((t1, t2))
        d = {t1: t1, t2: t2, t3: t3, "t3": t3}
        t3prime = jelly.unjelly(jelly.jelly(d))["t3"]
        self.assertIdentical(t3prime.other[0].other, t3prime.other[1].other)


    def test_classSecurity(self):
        """
        Test for class-level security of serialization.
        """
        taster = jelly.SecurityOptions()
        taster.allowInstancesOf(A, B)
        a = A()
        b = B()
        c = C()
        # add a little complexity to the data
        a.b = b
        a.c = c
        # and a backreference
        a.x = b
        b.c = c
        # first, a friendly insecure serialization
        friendly = jelly.jelly(a, taster)
        x = jelly.unjelly(friendly, taster)
        self.assertIsInstance(x.c, jelly.Unpersistable)
        # now, a malicious one
        mean = jelly.jelly(a)
        self.assertRaises(jelly.InsecureJelly, jelly.unjelly, mean, taster)
        self.assertIdentical(x.x, x.b, "Identity mismatch")
        # test class serialization
        friendly = jelly.jelly(A, taster)
        x = jelly.unjelly(friendly, taster)
        self.assertIdentical(x, A, "A came back: %s" % x)


    def test_unjellyable(self):
        """
        Test that if Unjellyable is used to deserialize a jellied object,
        state comes out right.
        """
        class JellyableTestClass(jelly.Jellyable):
            pass
        jelly.setUnjellyableForClass(JellyableTestClass, jelly.Unjellyable)
        input = JellyableTestClass()
        input.attribute = 'value'
        output = jelly.unjelly(jelly.jelly(input))
        self.assertEqual(output.attribute, 'value')
        self.assertIsInstance(output, jelly.Unjellyable)


    def test_persistentStorage(self):
        perst = [{}, 1]
        def persistentStore(obj, jel, perst = perst):
            perst[1] = perst[1] + 1
            perst[0][perst[1]] = obj
            return str(perst[1])

        def persistentLoad(pidstr, unj, perst = perst):
            pid = int(pidstr)
            return perst[0][pid]

        a = SimpleJellyTest(1, 2)
        b = SimpleJellyTest(3, 4)
        c = SimpleJellyTest(5, 6)

        a.b = b
        a.c = c
        c.b = b

        jel = jelly.jelly(a, persistentStore = persistentStore)
        x = jelly.unjelly(jel, persistentLoad = persistentLoad)

        self.assertIdentical(x.b, x.c.b)
        self.failUnless(perst[0], "persistentStore was not called.")
        self.assertIdentical(x.b, a.b, "Persistent storage identity failure.")


    def test_newStyleClassesAttributes(self):
        n = TestNode()
        n1 = TestNode(n)
        n11 = TestNode(n1)
        n2 = TestNode(n)
        # Jelly it
        jel = jelly.jelly(n)
        m = jelly.unjelly(jel)
        # Check that it has been restored ok
        self._check_newstyle(n, m)


    def _check_newstyle(self, a, b):
        self.assertEqual(a.id, b.id)
        self.assertEqual(a.classAttr, 4)
        self.assertEqual(b.classAttr, 4)
        self.assertEqual(len(a.children), len(b.children))
        for x, y in zip(a.children, b.children):
            self._check_newstyle(x, y)


    def test_referenceable(self):
        """
        A L{pb.Referenceable} instance jellies to a structure which unjellies to
        a L{pb.RemoteReference}.  The C{RemoteReference} has a I{luid} that
        matches up with the local object key in the L{pb.Broker} which sent the
        L{Referenceable}.
        """
        ref = pb.Referenceable()
        jellyBroker = pb.Broker()
        jellyBroker.makeConnection(StringTransport())
        j = jelly.jelly(ref, invoker=jellyBroker)

        unjellyBroker = pb.Broker()
        unjellyBroker.makeConnection(StringTransport())

        uj = jelly.unjelly(j, invoker=unjellyBroker)
        self.assertIn(uj.luid, jellyBroker.localObjects)



class ClassA(pb.Copyable, pb.RemoteCopy):
    def __init__(self):
        self.ref = ClassB(self)



class ClassB(pb.Copyable, pb.RemoteCopy):
    def __init__(self, ref):
        self.ref = ref



class CircularReferenceTestCase(unittest.TestCase):
    """
    Tests for circular references handling in the jelly/unjelly process.
    """

    def test_simpleCircle(self):
        jelly.setUnjellyableForClass(ClassA, ClassA)
        jelly.setUnjellyableForClass(ClassB, ClassB)
        a = jelly.unjelly(jelly.jelly(ClassA()))
        self.assertIdentical(a.ref.ref, a,
            "Identity not preserved in circular reference")


    def test_circleWithInvoker(self):
        class DummyInvokerClass:
            pass
        dummyInvoker = DummyInvokerClass()
        dummyInvoker.serializingPerspective = None
        a0 = ClassA()
        jelly.setUnjellyableForClass(ClassA, ClassA)
        jelly.setUnjellyableForClass(ClassB, ClassB)
        j = jelly.jelly(a0, invoker=dummyInvoker)
        a1 = jelly.unjelly(j)
        self.failUnlessIdentical(a1.ref.ref, a1,
            "Identity not preserved in circular reference")


    def test_set(self):
        """
        Check that a C{set} can contain a circular reference and be serialized
        and unserialized without losing the reference.
        """
        s = set()
        a = SimpleJellyTest(s, None)
        s.add(a)
        res = jelly.unjelly(jelly.jelly(a))
        self.assertIsInstance(res.x, set)
        self.assertEqual(list(res.x), [res])


    def test_frozenset(self):
        """
        Check that a C{frozenset} can contain a circular reference and be
        serializeserialized without losing the reference.
        """
        a = SimpleJellyTest(None, None)
        s = frozenset([a])
        a.x = s
        res = jelly.unjelly(jelly.jelly(a))
        self.assertIsInstance(res.x, frozenset)
        self.assertEqual(list(res.x), [res])

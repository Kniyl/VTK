# Copyright (c) Twisted Matrix Laboratories.
# See LICENSE for details.


from twisted.trial import unittest
import sys, os

from twisted.words.xish.domish import Element
from twisted.words.xish.xpath import XPathQuery
from twisted.words.xish import xpath

class XPathTest(unittest.TestCase):
    def setUp(self):
        # Build element:
        # <foo xmlns='testns' attrib1='value1' attrib3="user@host/resource">
        #     somecontent
        #     <bar>
        #        <foo>
        #         <gar>DEF</gar>
        #        </foo>
        #     </bar>
        #     somemorecontent
        #     <bar attrib2="value2">
        #        <bar>
        #          <foo/>
        #          <gar>ABC</gar>
        #        </bar>
        #     <bar/>
        #     <bar attrib4='value4' attrib5='value5'>
        #        <foo/>
        #        <gar>JKL</gar>
        #     </bar>
        #     <bar attrib4='value4' attrib5='value4'>
        #        <foo/>
        #        <gar>MNO</gar>
        #     </bar>
        #     <bar attrib4='value4' attrib5='value6'/>
        # </foo>
        self.e = Element(("testns", "foo"))
        self.e["attrib1"] = "value1"
        self.e["attrib3"] = "user@host/resource"
        self.e.addContent("somecontent")
        self.bar1 = self.e.addElement("bar")
        self.subfoo = self.bar1.addElement("foo")
        self.gar1 = self.subfoo.addElement("gar")
        self.gar1.addContent("DEF")
        self.e.addContent("somemorecontent")
        self.bar2 = self.e.addElement("bar")
        self.bar2["attrib2"] = "value2"
        self.bar3 = self.bar2.addElement("bar")
        self.subfoo2 = self.bar3.addElement("foo")
        self.gar2 = self.bar3.addElement("gar")
        self.gar2.addContent("ABC")
        self.bar4 = self.e.addElement("bar")
        self.bar5 = self.e.addElement("bar")
        self.bar5["attrib4"] = "value4"
        self.bar5["attrib5"] = "value5"
        self.subfoo3 = self.bar5.addElement("foo")
        self.gar3 = self.bar5.addElement("gar")
        self.gar3.addContent("JKL")
        self.bar6 = self.e.addElement("bar")
        self.bar6["attrib4"] = "value4"
        self.bar6["attrib5"] = "value4"
        self.subfoo4 = self.bar6.addElement("foo")
        self.gar4 = self.bar6.addElement("gar")
        self.gar4.addContent("MNO")
        self.bar7 = self.e.addElement("bar")
        self.bar7["attrib4"] = "value4"
        self.bar7["attrib5"] = "value6"

    def test_staticMethods(self):
        """
        Test basic operation of the static methods.
        """
        self.assertEqual(xpath.matches("/foo/bar", self.e),
                          True)
        self.assertEqual(xpath.queryForNodes("/foo/bar", self.e),
                          [self.bar1, self.bar2, self.bar4,
                           self.bar5, self.bar6, self.bar7])
        self.assertEqual(xpath.queryForString("/foo", self.e),
                          "somecontent")
        self.assertEqual(xpath.queryForStringList("/foo", self.e),
                          ["somecontent", "somemorecontent"])

    def test_locationFooBar(self):
        """
        Test matching foo with child bar.
        """
        xp = XPathQuery("/foo/bar")
        self.assertEqual(xp.matches(self.e), 1)

    def test_locationFooBarFoo(self):
        """
        Test finding foos at the second level.
        """
        xp = XPathQuery("/foo/bar/foo")
        self.assertEqual(xp.matches(self.e), 1)
        self.assertEqual(xp.queryForNodes(self.e), [self.subfoo,
                                                     self.subfoo3,
                                                     self.subfoo4])

    def test_locationNoBar3(self):
        """
        Test not finding bar3.
        """
        xp = XPathQuery("/foo/bar3")
        self.assertEqual(xp.matches(self.e), 0)

    def test_locationAllChilds(self):
        """
        Test finding childs of foo.
        """
        xp = XPathQuery("/foo/*")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar1, self.bar2,
                                                     self.bar4, self.bar5,
                                                     self.bar6, self.bar7])

    def test_attribute(self):
        """
        Test matching foo with attribute.
        """
        xp = XPathQuery("/foo[@attrib1]")
        self.assertEqual(xp.matches(self.e), True)

    def test_attributeWithValueAny(self):
        """
        Test find nodes with attribute having value.
        """
        xp = XPathQuery("/foo/*[@attrib2='value2']")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar2])

    def test_position(self):
        """
        Test finding element at position.
        """
        xp = XPathQuery("/foo/bar[2]")
        self.assertEqual(xp.matches(self.e), 1)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar1])

    test_position.todo = "XPath queries with position are not working."

    def test_namespaceFound(self):
        """
        Test matching node with namespace.
        """
        xp = XPathQuery("/foo[@xmlns='testns']/bar")
        self.assertEqual(xp.matches(self.e), 1)

    def test_namespaceNotFound(self):
        """
        Test not matching node with wrong namespace.
        """
        xp = XPathQuery("/foo[@xmlns='badns']/bar2")
        self.assertEqual(xp.matches(self.e), 0)

    def test_attributeWithValue(self):
        """
        Test matching node with attribute having value.
        """
        xp = XPathQuery("/foo[@attrib1='value1']")
        self.assertEqual(xp.matches(self.e), 1)

    def test_queryForString(self):
        """
        Test for queryForString and queryForStringList.
        """
        xp = XPathQuery("/foo")
        self.assertEqual(xp.queryForString(self.e), "somecontent")
        self.assertEqual(xp.queryForStringList(self.e),
                          ["somecontent", "somemorecontent"])

    def test_queryForNodes(self):
        """
        Test finding nodes.
        """
        xp = XPathQuery("/foo/bar")
        self.assertEqual(xp.queryForNodes(self.e), [self.bar1, self.bar2,
                                                     self.bar4, self.bar5,
                                                     self.bar6, self.bar7])

    def test_textCondition(self):
        """
        Test matching a node with given text.
        """
        xp = XPathQuery("/foo[text() = 'somecontent']")
        self.assertEqual(xp.matches(self.e), True)

    def test_textNotOperator(self):
        """
        Test for not operator.
        """
        xp = XPathQuery("/foo[not(@nosuchattrib)]")
        self.assertEqual(xp.matches(self.e), True)

    def test_anyLocationAndText(self):
        """
        Test finding any nodes named gar and getting their text contents.
        """
        xp = XPathQuery("//gar")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.gar1, self.gar2,
                                                     self.gar3, self.gar4])
        self.assertEqual(xp.queryForStringList(self.e), ["DEF", "ABC",
                                                          "JKL", "MNO"])

    def test_anyLocation(self):
        """
        Test finding any nodes named bar.
        """
        xp = XPathQuery("//bar")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar1, self.bar2,
                                                     self.bar3, self.bar4,
                                                     self.bar5, self.bar6,
                                                     self.bar7])

    def test_anyLocationQueryForString(self):
        """
        L{XPathQuery.queryForString} should raise a L{NotImplementedError}
        for any location.
        """
        xp = XPathQuery("//bar")
        self.assertRaises(NotImplementedError, xp.queryForString, None)

    def test_andOperator(self):
        """
        Test boolean and operator in condition.
        """
        xp = XPathQuery("//bar[@attrib4='value4' and @attrib5='value5']")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar5])

    def test_orOperator(self):
        """
        Test boolean or operator in condition.
        """
        xp = XPathQuery("//bar[@attrib5='value4' or @attrib5='value5']")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar5, self.bar6])

    def test_booleanOperatorsParens(self):
        """
        Test multiple boolean operators in condition with parens.
        """
        xp = XPathQuery("""//bar[@attrib4='value4' and
                                 (@attrib5='value4' or @attrib5='value6')]""")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar6, self.bar7])

    def test_booleanOperatorsNoParens(self):
        """
        Test multiple boolean operators in condition without parens.
        """
        xp = XPathQuery("""//bar[@attrib5='value4' or
                                 @attrib5='value5' or
                                 @attrib5='value6']""")
        self.assertEqual(xp.matches(self.e), True)
        self.assertEqual(xp.queryForNodes(self.e), [self.bar5, self.bar6, self.bar7])

FluDAG Interface Testing
===========================

These tests provide the ability to test individual geometry queries of
the FluDAG interface.  They will be divided into sections, each
section representing one high level query type, with individual tests
for the different use cases and edge cases.

Distance to Boundary
--------------------

This is the most fundamental geoemtry query.  This also have the most
use and edge cases in the DAGMC paradigm due to the use of the ray
history (aka "tail") to provide robustness and efficiency.

Point Inclusion & Zone Search
--------------------------------

The DAGMC interface historically provides only the means to test
whether or not a point is inside a given zone.  FluDAG requires the
ability to identify which zone houses a point, so a search through
zones will also be necessary.

Boundary Crossings: Normals and Next Zone
--------------------------------------------

Utility functions require the normal to a surface at a point and
information about what zone is entered upon crossing a surface.

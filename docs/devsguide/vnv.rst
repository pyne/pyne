Verification and Validation
===========================

A PyNE module can be declared Verification and Validation (V&V) compliant provided 

Once a PyNE module is verified and validated 

American Society of Mechanical Engineers (ASME) NQA-1-2008/NQA-1a-2009 Part II Subpart 2.7 Requirement 302 "Otherwise Acquired Software"



Criteria for declaring and module V&V compliant
------------------------------------------------

1. The module has a theory manual entry that describes or cites all non-trivial mathematics or physics that occur within a module, and also lists the assumptions and limitations of the module.
2. The module is fully unit tested, with test coverage that spans the entire set of capabilities described in the theory manual.  If the expected results of unit tests require non-trivial calculations, sample calculations should be provided in the theory manual. In other words, any numbers that appear in the unit tests must be reproducible.
3. The module has a user manual entry that describes how to use its functionality within the design basis.
4. If modules contain significant physics (e.g. transport, transmutation), the module must be benchmarked against experimental results, with appropriate uncertainty and sensitivity analysis. In addition, regression tests must be incorporated into the CI.
5. When removal of an import guard is proposed in a pull-request, both the requester and reviewer should be prepared to justify the decision to merge in the event that the module is subject to the scrutiny of a third party.

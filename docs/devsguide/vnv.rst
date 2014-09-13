Verification and Validation
===========================

Software Quality Assurance (QA) standards describe requirements for both software
products (code, tests, documentation) and also software development practices.
The software development practices within PyNE can be mapped to those
prescribed by industry standards [1]. In regard to software products, due to PyNE's status as an open-source,
community-developed and ongoing project, the degree to which different areas of
the codebase are rigorously verified, validated, and documented is not expected to be
homogeneous. For this reason, PyNE code is declared Verification and
Validation (V&V) complaint on a per-module basis. Modules that are not V&V
compliant must contain Python import warnings.

In order to be declared V&V compliant, a module must meet the criteria
described below. These V&V criteria are internal to PyNE, but they were
developed by adapting other QA criteria important within the nuclear
engineering field. The Nuclear Regulatory Commission (NRC) endorses the
American Society of Mechanical Engineers (ASME) NQA-1-2008/NQA-1a-2009 quality
assurance standard for the licensing of power plants and reprocessing
facilities. Since PyNE, for the most part, is a library intended to be used as
a component within an end-user application, those who seek to use PyNE code for
NQA-1 compliant purposes will likely evoke Part II
Subpart 2.7 Requirement 302 "Otherwise Acquired Software" of NQA-1-2008. The PyNE V&V
criteria ensure that Requirement 302 criteria are addressed and readily
available within PyNE documentation. In addition, the PyNE V&V criteria are
intended to correspond roughly with the Quality Rigor Level 1 Verification and
Validation requirements as defined by the Nuclear Energy Modeling and Simulation
(NEAMS) Program toolkit verification and validation plan `[pdf]
<http://www.energy.gov/sites/prod/files/2013/09/f2/NEAMS%20Software%20Verification%20and%20Validation%20Plan%20Requirements%20Version%200.pdf>`_.
This will facilitate the process of formally declaring PyNE modules compliant
with these standards, if desired by some party.

Criteria for declaring a module V&V compliant
------------------------------------------------

1. The module has a theory manual entry that describes or cites all non-trivial mathematics and/or physics that occur within a module, and also lists the assumptions and limitations of the module.
2. The module is fully unit tested, with test coverage that spans the entire set of capabilities described in the theory manual.  If the expected results of unit tests require non-trivial calculations, sample calculations should be provided in the theory manual. In other words, any numbers that appear in the unit tests must be reproducible.
3. The module has a user manual entry that describes how to use its functionality within the design basis.
4. If modules contain significant physics (e.g. transport, transmutation), the module must be benchmarked against experimental results, with appropriate uncertainty and sensitivity analysis. In addition, regression tests must be incorporated into the CI.
5. When removal of an import guard is proposed in a pull-request, both the requester and reviewer should be prepared to justify the decision to merge, in the event that the module falls under scrutiny.

References
-----------
[1] E. Biondo, A. Scopatz, M. Gidden, R. Slaybaugh, C. Bates, P. P.H. Wilson, *Quality Assurance within the PyNE Nuclear Toolkit*, Transactions of the American Nuclear Society, 111, *In Press*. 


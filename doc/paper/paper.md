---
title: 'Pencil Code, a modular MPI code for partial differential equations: multipurpose and multiuser-maintained'
tags:
 - Fortran90
 - fluid dynamics
 - magnetohydrodynamics
 - Python
 - IDL
 - astrophysics
 - radiation
 - inertial particles
 - combustion
authors:
 - name: The Pencil Code Collaboration
   affiliation: 1
 - name: Axel Brandenburg
   affiliation: "1, 2, 3"
   orcid: 0000-0002-7304-021X
 - name: Anders Johansen
   affiliation: 4
   orcid: 0000-0002-5893-6165
 - name: Philippe A. Bourdin
   affiliation: 5
   orcid: 0000-0002-6793-601X
 - name: Wolfgang Dobler
   affiliation: 6
 - name: Wladimir Lyra
   affiliation: 7
   orcid: 0000-0002-3768-7542
 - name: Matthias Rheinhardt
   affiliation: 8
 - name: Nils Erland L. Haugen
   affiliation: "9, 10, 1"
   orcid: 0000-0002-9184-8722
 - name: Antony Mee
   affiliation: 11
 - name: Frederick Gent
   affiliation: "8, 12"
   orcid: 0000-0002-1331-2260
 - name: Natalia Babkovskaia
   affiliation: 13
 - name: Chao-Chin Yang
   affiliation: 14
   orcid: 0000-0003-2589-5034
 - name: Simon Candelaresi
   affiliation: 15
   orcid: 0000-0002-7666-8504
 - name: Jörn Warnecke
   affiliation: 16
   orcid: 0000-0002-9292-4600
 - name: Petri Käpylä
   affiliation: 17
   orcid: 0000-0001-9619-0053
 - name: Piyali Chatterjee
   affiliation: 18
   orcid: 0000-0002-0181-2495
 - name: Xiang-Yu Li
   affiliation: 19
   orcid: 0000-0002-5722-0018 
 - name: Jonas Krüger
   affiliation: "9, 10"
 - name: Jørgen R. Aarnes
   affiliation: 10
   orcid: 0000-0002-5899-2597
 - name: Graeme Sarson
   affiliation: 12
 - name: Jennifer Schober
   affiliation: 20
   orcid: 0000-0001-7888-6671
 - name: Raphaël Plasson
   affiliation: 21
   orcid: 0000-0003-2319-1463
 - name: Christer Sandin
   affiliation: 1
   orcid: 0000-0002-6370-5505
 - name: Luiz Felippe S Rodrigues
   affiliation: "12, 22"
   orcid: 0000-0002-3860-0525
 - name: Add Yourself
   affiliation: 23

affiliations:
 - name: Nordita, KTH Royal Institute of Technology and Stockholm University
   index: 1
 - name: Department of Astronomy, Stockholm University
   index: 2
 - name: McWilliams Center for Cosmology & Department of Physics, Carnegie Mellon University
   index: 3
 - name: GLOBE Institute, University of Copenhagen
   index: 4
 - name: Space Research Institute, Graz
   index: 5
 - name: Bruker, Potsdam
   index: 6
 - name: New Mexico State University, Department of Astronomy, PO Box 30001, MSC 4500 Las Cruces, NM 88001
   index: 7
 - name: Astroinformatics, Department of Computer Science, Aalto University
   index: 8
 - name: SINTEF Energy Research, Trondheim, Norway
   index: 9
 - name: Norwegian University of Science and Technology
   index: 10
 - name: Barclays, London
   index: 11
 - name: School of Mathematics, Statistics and Physics, Newcastle University 
   index: 12
 - name: No current affiliation
   index: 13
 - name: University of Nevada, Las Vegas
   index: 14
 - name: School of Mathematics and Statistics, University of Glasgow
   index: 15
 - name: Max Planck Institute for Solar System Research 
   index: 16
 - name: Institute for Astrophysics, University of Göttingen
   index: 17
 - name: Indian Institute of Astrophysics, Bengaluru-560034, India
   index: 18
 - name: Pacific Northwest National Laboratory
   index: 19
 - name: Laboratoire d'Astrophysique, EPFL, CH-1290 Sauverny, Switzerland
   index: 20
 - name: Avignon Université, France
   index: 21
 - name: Radboud University, Netherlands
   index: 22
 - name: Currently The Last One
   index: 23

date: 20 August 2020
bibliography: paper.bib
---

# Summary

The Pencil Code is a highly modular physics-oriented simulation code
that can be adapted to a wide range of applications.
It is primarily designed to solve partial differential equations (PDEs)
of compressible hydrodynamics and has lots of add-ons ranging from
astrophysical magnetohydrodynamics (MHD) [@2010ascl.soft10060B] to
meteorological cloud microphysics [@2017JAMES.9.1116L] and engineering
applications in combustion [@2011JCoPh.230.1B].
Nevertheless, the framework is general and can also be applied to
situations not related to hydrodynamics or even PDEs, for example when
just the message passing interface or input/output strategies of the
code are to be used.
The code can also solve for Lagrangian (inertial and noninertial)
particles, their coagulation and condensation, as well as their
interaction with the fluid.
A related module has also been adapted to perform ray tracing
to solve the eikonal equation.

The code is being used for Cartesian, cylindrical, and spherical geometries,
but further extensions are possible.
One can choose between different time stepping schemes and different
spatial derivative operators.
High-order first and second derivatives are used to deal with weakly 
compressible turbulent flows.
There are also different diffusion operators to allow for both direct numerical
simulations (DNS) and various types of large-eddy simulations (LES).

# High-level functionality

An idea about the range of available modules can be obtained by inspecting
the examples under pencil-code/samples/.
Those are low resolution versions related to applications published in the literature.
Some of the run directories of actual production runs are published through Zenodo.
Below a list of method papers that describe the various applications and tests:

* Coagulation and condensation in turbulence [@2008A&A.486.597J; @2017JAMES.9.1116L]
* Radiative transfer [@2006A&A.448.731H; @2014A&A.571A.68B; @2020GApFD.114.162B]
* Chiral magnetic effect in relativistic plasmas [@2018ApJ.858.124S]
* Primordial gravitational waves [@2020GApFD.114.130R]
* Modeling homochirality at the origin of life [@2004IJAsB.3.209B; @2019OLEB.49.49B]
* Modeling of patterned photochemical systems [@2012ChemEurJ]
* Gaseous combustion and detonation [@2011JCoPh.230.1B; @Zhang_etal_2020comb]
* Burning particles, resolved or unresolved [@2020GApFD.114.58Q; @2017CNF.185a160]
* Flows around immersed solid objects [@2019IJCFD.33.43A; @2020GApFD.114.35A; @2010JFM.661a239]
* Test-field method for turbulent MHD transport [@2010A&A.520A.28R; @2010PhST.142a4028B; @2018A&A.609A.51W]
* Spherical shell dynamos and convection [@2009ApJ.697.923M; @2020GApFD.114.8K]
* Boris correction for coronal physics [@2020GApFD.114.213C]
* Thermal instability and mixing [@2012ApJ.758.48Y]
* Implicit solver for temperature [@2008A&A.484.29G]
* Dust-gas dynamics with mutual drag interaction [@2007ApJ.662.613Y; @2016ApJS.224.39Y]

# Statement of need and purpose of software

The code provides an easily adaptable tool for solving both standard
MHD equations as well as others, such as the test-field equations.
Significant amounts of runtime diagnostics 
as well as Python and IDL libraries for post-processing are available.

Among the currently 83 developers with check-in permission, there are
currently 16 owners who can give others check-in permission.
Of the developers, 34 have done more than 34 commits.
There are others with fewer commits who have contributed with more than
5000 lines to the code and also contributed significantly to the code.
Users have access to the latest development version and can ask to
join the circle of developers by contacting one of the owners.

Every revision on GitHub is verified on 9 tests on travis-ci.com.
The current version is also automatically being tested on 59 hourly
tests and on 79 daily tests.
Continuous progress on the code is driven by the research of
individual developers.

Further developments and interactions between developers and users are
being promoted through annual user meetings since 2004 and a newsletters
since 2020.
Since 2016, a steering committee of five elected owners reviews the
progress and can take decisions of general concern to the Pencil Code
community.

# Ongoing research using the Pencil Code

Current research topics with the code includes
interstellar and intercluster medium as well as early Universe;
small-scale dynamos and reconnection;
primordial magnetic fields and decaying turbulence;
gravitational waves from turbulent sources;
planet formation and inertial particles;
accretion discs and shear flows;
coronal heating and coronal mass ejections;
helical dynamos, helical turbulence, and catastrophic quenching;
helioseismology;
strongly stratified MHD turbulence and negative effective magnetic pressure instability;
convection in Cartesian domains;
global convection and dynamo simulations;
turbulent transport and test-field methods;
hydrodynamic and MHD instabilities;
chiral MHD;
hydrodynamic and MHD turbulence;
turbulent gaseous and solid combustion, particle clustering and deposition on solid walls,
front propagation, radiation & ionization.
As of July 2020, 564 papers have been published that acknowledge use of
the Pencil Code [@zenodo.3466444].

# Key references

The Pencil Code is unique in two ways:
the high level of flexibility and modularity, and the way it is organized
(open source, distributed ownership, openness of development version).

Other software addressing related needs include: 
Athena, CO5BOLD, ENZO, MuRAM, NIRVANA, Stagger, ZEUS, and several other LES codes.
There are also several other engineering DNS codes such as
the Sandia-3-Dimensional (S3D) high-order compressible
optimized for combustion, which is not open source, however.
In addition, there is the Dedalus code, which allows one to
program the equations in symbolic form.

Some recent research areas that made use of the Pencil Code, as
evidenced by the aforementioned document listing all those papers
[@zenodo.3466444], include:

* Flows around immersed solid objects [@2010JFM.661a239]
* Particle clustering in supersonic and subsonic turbulence [@2019MNRAS.483.5623M; @Karchniwy_etal_2019],
* Cloud microphysics [@2017JAMES.9.1116L],
* Planet and planetesimal formation [@2007Natur.448.1022J; @2007ApJ.670.805O; @2009A&A.497.869L],
* Global simulations of debris disks [@2013Natur.499.184L],
* Stratified shearing box simulations, also with dust [@2011ApJ.740.18O; @2018ApJ.861.47S; @2018ApJ.868.27Y].
* Supernova-driven turbulence [@2013MNRAS.432.1396G],
* Solar dynamo and sunspots [@2005ApJ.625.539B; @2007ApJ.669.1390H],
* Solar corona above active regions [@2011A&A.530A.112B; @2013A&A.555A.123B; @2016PhRvL.116j1101C],
* Fully convective star in a box [@2006ApJ.638.336D],
* Dynamo wave in spherical shell convection [@2012ApJ.755L.22K; @2014ApJ.796L.12W],
* MHD turbulence and cascades [@2004PhRvE.70a6308H],
* Turbulent diffusivity quenching with test fields [@2008ApJ.676.740B; @2014ApJ.795.16K],

# Acknowledgements

We acknowledge contributions from all submitters and their supporting
funding agencies.
In particular, we mention the Swedish Research Council,
grants 2012-5797, 2013-03992, 2017-03865, and 2019-04234,
the National Science Foundation under the grant AAG-1615100,
the FRINATEK grant 231444 under the Research Council of Norway, SeRC,
the grant "Bottlenecks for particle growth in turbulent aerosols"
from the Knut and Alice Wallenberg Foundation, Dnr.\ KAW 2014.0048,
and the University of Colorado through its support of the
George Ellery Hale visiting faculty appointment,

# References


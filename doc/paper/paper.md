---
title: 'The Pencil Code, a modular MPI code for partial differential equations and particles: multipurpose and multiuser-maintained'

# The Pencil Code is used and developed by the 37 authors, who define the Pencil Code Collaboration.
# About half of the currently 560 papers that acknowledge the code are by others who picked up the code, which has been public since 2001.

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
 - name: Sven Bingert
   affiliation: 9
   orcid: 0000-0001-9547-1582
 - name: Nils Erland L. Haugen
   affiliation: "10, 11, 1"
   orcid: 0000-0002-9184-8722
 - name: Antony Mee
   affiliation: 12
 - name: Frederick Gent
   affiliation: "8, 13"
   orcid: 0000-0002-1331-2260
 - name: Natalia Babkovskaia
   affiliation: 14
 - name: Chao-Chin Yang
   affiliation: 15
   orcid: 0000-0003-2589-5034
 - name: Tobias Heinemann
   affiliation: 16
   orcid: 0000-0002-2476-9733
 - name: Boris Dintrans
   affiliation: 17
 - name: Dhrubaditya Mitra
   affiliation: 1
   orcid: 0000-0003-4861-8152
 - name: Simon Candelaresi
   affiliation: 18
   orcid: 0000-0002-7666-8504
 - name: Jörn Warnecke
   affiliation: 19
   orcid: 0000-0002-9292-4600
 - name: Petri J. Käpylä
   affiliation: 20
   orcid: 0000-0001-9619-0053
 - name: Andreas Schreiber
   affiliation: 14
 - name: Piyali Chatterjee
   affiliation: 21
   orcid: 0000-0002-0181-2495
 - name: Maarit J. Käpylä
   affiliation: "8, 19"
   orcid: 0000-0002-9614-2200
 - name: Xiang-Yu Li
   affiliation: 1
   orcid: 0000-0002-5722-0018 
 - name: Jonas Krüger
   affiliation: "10, 11"
   orcid: 0000-0001-8036-0695
 - name: Jørgen R. Aarnes
   affiliation: 11
   orcid: 0000-0002-5899-2597
 - name: Graeme R. Sarson
   affiliation: 13
   orcid: 0000-0001-6774-9372
 - name: Jeffrey S. Oishi
   affiliation: 22
   orcid: 0000-0001-8531-6570
 - name: Jennifer Schober
   affiliation: 23
   orcid: 0000-0001-7888-6671
 - name: Raphaël Plasson
   affiliation: 24
   orcid: 0000-0003-2319-1463
 - name: Christer Sandin
   affiliation: 1
   orcid: 0000-0002-6370-5505
 - name: Ewa Karchniwy
   affiliation: "11, 25"
   orcid: 0000-0001-6709-1160
 - name: Luiz Felippe S. Rodrigues
   affiliation: "13, 26"
   orcid: 0000-0002-3860-0525
 - name: Alexander Hubbard
   affiliation: 27
 - name: Gustavo Guerrero
   affiliation: 28
   orcid: 0000-0002-2671-8796
 - name: Andrew Snodin
   affiliation: 13
 - name: Illa R. Losada
   affiliation: 1
   orcid: 0000-0002-0416-7516
 - name: Johannes Pekkilä
   affiliation: 8
   orcid: 0000-0002-1974-7150
 - name: Chengeng Qian
   affiliation: 29
   orcid: 0000-0002-5560-5475

affiliations:
 - name: Nordita, KTH Royal Institute of Technology and Stockholm University, Sweden
   index: 1
 - name: Department of Astronomy, Stockholm University, Sweden
   index: 2
 - name: McWilliams Center for Cosmology & Department of Physics, Carnegie Mellon University, PA, USA
   index: 3
 - name: GLOBE Institute, University of Copenhagen, Denmark
   index: 4
 - name: Space Research Institute, Graz, Austria
   index: 5
 - name: Bruker, Potsdam, Germany
   index: 6
 - name: New Mexico State University, Department of Astronomy, Las Cruces, NM, USA
   index: 7
 - name: Astroinformatics, Department of Computer Science, Aalto University, Finland
   index: 8
 - name: Gesellschaft für wissenschaftliche Datenverarbeitung mbH Göttingen, Germany
   index: 9
 - name: SINTEF Energy Research, Trondheim, Norway
   index: 10
 - name: Norwegian University of Science and Technology, Norway
   index: 11
 - name: Bank of America Merrill Lynch, London, UK
   index: 12
 - name: School of Mathematics, Statistics and Physics, Newcastle University, UK
   index: 13
 - name: No current affiliation
   index: 14
 - name: University of Nevada, Las Vegas, USA
   index: 15
 - name: Niels Bohr International Academy, Denmark
   index: 16
 - name: CINES, Montpellier, France
   index: 17
 - name: School of Mathematics and Statistics, University of Glasgow, UK
   index: 18
 - name: Max Planck Institute for Solar System Research, Germany
   index: 19
 - name: Institute for Astrophysics, University of Göttinge, Germany
   index: 20
 - name: Indian Institute of Astrophysics, Bengaluru, India
   index: 21
 - name: Department of Physics & Astronomy, Bates College, ME, USA
   index: 22
 - name: Laboratoire d'Astrophysique, EPFL, Sauverny, Switzerland
   index: 23
 - name: Avignon Université, France
   index: 24
 - name: Institute of Thermal Technology, Silesian University of Technology, Poland
   index: 25
 - name: Radboud University, Netherlands
   index: 26
 - name: Department of Astrophysics, American Museum of Natural History, NY, USA
   index: 27
 - name: Physics Department, Universidade Federal de Minas Gerais, Belo Horizonte, Brazil
   index: 28
 - name: State Key Laboratory of Explosion Science and Technology, Beijing Institute of Technology, China
   index: 29

date: 17 September 2020
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
The code can also evolve Lagrangian (inertial and noninertial)
particles, their coagulation and condensation, as well as their
interaction with the fluid.
A related module has also been adapted to perform ray tracing
and to solve the eikonal equation.

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

* Coagulation and condensation in turbulence [@2008A&A.486.597J; @2017JAMES.9.1116L],
* Radiative transfer [@2006A&A.448.731H; @2014A&A.571A.68B; @2020GApFD.114.162B],
* Chiral magnetic effect in relativistic plasmas [@2018ApJ.858.124S],
* Primordial gravitational waves [@2020GApFD.114.130R],
* Modeling homochirality at the origin of life [@2004IJAsB.3.209B; @2019OLEB.49.49B],
* Modeling of patterned photochemical systems [@2012ChemEurJ],
* Gaseous combustion and detonation [@2011JCoPh.230.1B; @Zhang_etal_2020comb; @2017CNF.185a160],
* Burning particles, resolved or unresolved [@2020GApFD.114.58Q],
* Flows around immersed solid objects [@2019IJCFD.33.43A; @2020GApFD.114.35A; @2010JFM.661a239],
* Test-field method for turbulent MHD transport [@2010A&A.520A.28R; @2010PhST.142a4028B; @2018A&A.609A.51W],
* Mean-field MHD [@2013SoPh.287.293K; @2013A&A.556A.106J],
* Spherical shell dynamos and convection [@2009ApJ.697.923M; @2020GApFD.114.8K],
* Boris correction for coronal physics [@2020GApFD.114.213C],
* Thermal instability and mixing [@2012ApJ.758.48Y],
* Implicit solver for temperature [@2008A&A.484.29G],
* Dust-gas dynamics with mutual drag interaction [@2007ApJ.662.613Y; @2016ApJS.224.39Y].

# Statement of need and purpose of software

The code is an easily adaptable tool for solving both standard
MHD equations as well as others, such as the test-field equations.
Significant amounts of runtime diagnostics 
as well as Python and IDL libraries for post-processing are available.

Among the currently 83 developers with check-in permission, there are
currently 18 owners who can give others check-in permission.
Of the developers, 34 have done more than 34 commits.
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

Current research includes topics from stellar physics, interstellar and intercluster medium, the early universe,
as well as from meteorology and engineering:
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
hydrodynamic and MHD instabilities and turbulence;
chiral MHD;
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
Sandia-3-Dimensional (S3D), a high-order compressible code,
optimized for combustion, which is not open source, however.
In addition, there are frameworks like Dedalus or Cactus,
which allow one to program the equations in symbolic form.

Some recent research areas that made use of the Pencil Code, as
evidenced by the aforementioned document listing all those papers
[@zenodo.3466444], include:

* Flows around immersed solid objects [@2010JFM.661a239],
* Particle clustering in supersonic and subsonic turbulence [@2019MNRAS.483.5623M; @Karchniwy_etal_2019],
* Cloud microphysics [@2017JAMES.9.1116L],
* Planet and planetesimal formation [@2007Natur.448.1022J; @2007ApJ.670.805O; @2009A&A.497.869L],
* Global simulations of debris disks [@2013Natur.499.184L],
* Stratified shearing box simulations, also with dust [@2011ApJ.740.18O; @2018ApJ.861.47S; @2018ApJ.868.27Y],
* Supernova-driven turbulence [@2013MNRAS.432.1396G],
* Solar dynamo and sunspots [@2005ApJ.625.539B; @2007ApJ.669.1390H],
* Solar corona above active regions [@2011A&A.530A.112B; @2013A&A.555A.123B; @2016PhRvL.116j1101C],
* Fully convective star in a box [@2006ApJ.638.336D],
* Dynamo wave in spherical shell convection [@2012ApJ.755L.22K; @2014ApJ.796L.12W],
* Convection with Kramers opacity law [@2017ApJ.845.23K; @2019A&A.631.122K; @2020GApFD.114.8K],
* MHD turbulence and cascades [@2004PhRvE.70a6308H],
* Turbulent diffusivity quenching with test fields [@2008ApJ.676.740B; @2014ApJ.795.16K].

# Acknowledgements

We acknowledge contributions from all submitters and their supporting
funding agencies.
In particular, we mention the ERC Advanced Grant on Astrophysical Dynamos
(No 227952), the Swedish Research Council,
grants 2012-5797, 2013-03992, 2017-03865, and 2019-04234,
the National Science Foundation under the grant AAG-1615100,
the FRINATEK grant 231444 under the Research Council of Norway, SeRC,
the grant "Bottlenecks for particle growth in turbulent aerosols"
from the Knut and Alice Wallenberg Foundation, Dnr.\ KAW 2014.0048,
the ReSoLVE Centre of Excellence (grant number 307411),
the research project ‘Gaspro’, financed by the Research Council of
Norway (267916),
the European Research Council (ERC) under the European Union's
Horizon 2020 research and innovation programme (Project UniSDyn,
grant agreement n:o 818665), and the Deutsche Forschungsgemeinschaft
(DFG) Heisenberg programme grant KA 4825/2-1.

# References


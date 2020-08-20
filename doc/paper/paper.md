---
title: 'Pencil Code, a modular MPI code for partial differential equation: multipurpose and multiuser-maintained'
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
 - name: Philippe A. Bourdin
   affiliation: 5
 - name: Wolfgang Dobler
   affiliation: 6
 - name: Wladimir Lyra
   affiliation: 7
 - name: Add Yourself
   affiliation: 8
affiliations:
 - name: Nordita, KTH Royal Institute of Technology and Stockholm University
   index: 1
 - name: Department of Astronomy, Stockholm University
   index: 2
 - name: McWilliams Center for Cosmology & Department of Physics, Carnegie Mellon University
   index: 3
 - name: New Mexico State Univ.
   index: 4
 - name: Four
   index: 5
 - name: Five
   index: 6
 - name: Six
   index: 7
 - name: Seven
   index: 8
 - name: Eight
date: 20 August 2020
bibliography: paper.bib
---

# Summary

The Pencil Code is highly modular and can be adapted to a wide
range of applications.
It is primarily designed to solving partial differential equations
(PDEs) of compressible magnetohydrodynamics [@2010ascl.soft10060B],
but the framework is general and can also be applied in situations not
related to PDEs, for example when just the message passing interface or
input/output strategies of the code are to be used.

The code is being used for Cartesian, cylindrical, and spherical geometries,
but further extensions are possible.
A preliminary implementation of a Yin--Yang mesh is also in place.
One can choose between different time stepping schemes and different
spatial derivative operators.
There are also different diffusion operators to allow for both direct numerical
simulations (DNS) and various types of large-eddy simulations (LES).

# High-level functionality

Coagulation and condensation in turbulence [@2017JAMES...9.1116L],
radiative transfer [@2006A&A...448..731H,@2020GApFD.114..162B],
chiral magnetic effect in relativistic plasmas [@2018ApJ...858..124S],
primordial gravitational waves [@2020GApFD.114..130R],
combustion and detonation [@2011JCoPh.230....1B,@Zhang_etal_2020comb,@2020GApFD.114...58Q],
test-field method for turbulent MHD transport coefficients [@2010PhST..142a4028B],
spherical shell convection [@2020GApFD.114....8K].

* Coagulation and condensation in turbulence [@2017JAMES...9.1116L]
* Radiative transfer [@2006A&A...448..731H,@2020GApFD.114..162B]
* Chiral magnetic effect in relativistic plasmas [@2018ApJ...858..124S]
* Primordial gravitational waves [@2020GApFD.114..130R]
* Combustion and detonation [@2011JCoPh.230....1B,@Zhang_etal_2020comb,@2020GApFD.114...58Q]
* Test-field method for turbulent MHD transport coefficients [@2010PhST..142a4028B]
* Spherical shell convection [@2020GApFD.114....8K]

# Statement of need and purpose of software

The code provides an easily adaptable tool for solving both well
established equations and new ones, such as the test-field equations.
Significant amounts of runtime diagnostics is available.
Users have access to the latest development version and can
join the circle of developers.
Every revision on GitHub is verified on 12 tests on travis-ci.com.
The current version is also automatically being tested on 59 hourly
tests and on 79 daily tests.
Continuous progress on the code is driven by the research of the
individual developers.
Among the currently 86 developers with check-in permission, there are
currently 16 owners who can give others check-in permission.
Of the developers, 34 have done more than 34 comments.
Further developments and interactions between developers and users is
being promoted through annual user meetings since 2004 and a newsletters
since 2020.
Since 2016, a steering committee of five elected owners reviews the
progress and can take decisions of general concern to the Pencil Code
community.

# Ongoing research using the Pencil Code

Current research with the code includes solar and stellar convection,
planetesimal formation, coagulation and condensation in raindrop and dust formation,
turbulent combustion, and gravitational waves from the early universe.
As of July 2020, 564 papers have been published that acknowledge use of
the Pencil Code [@zenodo.3466444].

# Key references

* Planet formation [@2007Natur.448.1022J]
* Solar dynamo [@2005ApJ...625..539B]
* MHD turbulence and cascades [@2004PhRvE..70a6308H]
* Fully convective star in a box [@2006ApJ...638..336D]
* Dynamo wave in spherical shell convection [@2012ApJ...755L..22K]
* Turbulent diffusivity quenching with test fields [@2008ApJ...676..740B,@2014ApJ...795...16K]
* Global simulations of debris disks [@2013Natur.499..184L]
* Supernova-driven turbulence [@2013MNRAS.432.1396G]
* Clustering in in supersonic turbulent molecular clouds [@2019MNRAS.483.5623M]
* Solar corona above active regions [@2013A&A...555A.123B]
* Stratified shearing box simulations [@2011ApJ...740...18O]

# Acknowledgements

We acknowledge contributions from all submitters and their supporting
funding agencies.
This work was supported in part through the Swedish Research Council,
grants 2012-5797, 2013-03992, 2017-03865, and 2019-04234,
the National Science Foundation under the grant AAG-1615100,
the FRINATEK grant 231444 under the Research Council of Norway, SeRC,
the grant ``Bottlenecks for particle growth in turbulent aerosols''
from the Knut and Alice Wallenberg Foundation, Dnr.\ KAW 2014.0048,
and the University of Colorado through its support of the
George Ellery Hale visiting faculty appointment,

# References


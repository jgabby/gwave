# Introduction

The FCC's [OET Report 86](oet_r86-1.pdf) outlined a computer program to automate the generation of Groundwave propagation charts.  Code was provided in Fortran.  AM program methodology was also included in International Agreements with [Canada](https://www.fcc.gov/broadcast-agreements-canada), [Mexico](https://www.fcc.gov/broadcast-agreements-mexico) and [Region 2](rj81.pdf)

Subsequent internal implementations of the code were automatically converted from Fortran to C but included non-public libraries.  The program also makes extensive use of Complex numbers, which are not neccesarily included in some current programming languages.

# Purpose
This project reviewed the historical documentation and existing code base to develop a javascript module capable of generating groundwave charts.  These charts, in the form of JSON tables, can be used in AM groundwave propagation studies.

To handle complex numbers, custom functions have been included.

# Installation

Clone the repository, cd into the folder then run:
```
npm install
node server.js
```
A message will appear showing that the service is available on port 3000.


# Usage
The only input for this function is the frequency to study in KHz:
Example: http://localhost:3000/?frequency=550


# Reflections
It's worth noting that while the gwave program is capable of directly solving for any specific frequency, distance, and conductivity, the FCC chose to generate tables and then interpolate within the tables.  
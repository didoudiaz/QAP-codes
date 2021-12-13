# QAP-codes
This project provides several codes for the Quadratic Assignment Problems (solvers, utilities,...).

Includes an extended Extremal Optimization (EO) procedure for QAP.

## Usage

Should compile on linux simply with:

        cd src
        make

to see help use -h (generally works)

        eo-qap -h
        rots-qap -h

then run, e.g. with default parameter values and some verbosity:

        eo-qap -v 1 ~/QAP-instances/QAPLIB/tai40a.qap 
        eo-rots -v 1 ~/QAP-instances/QAPLIB/tai40a.qap 

## See also

[codes on QAPLIB](http://www.mgi.polymtl.ca/anjos/qaplib/codes.html)

[A large collection of QAP instances](https://github.com/didoudiaz/QAP-instances)

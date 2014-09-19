amsstats
=========

Statistics tools for determining Anisotropy of Magnetic Susceptibility

Installation
------------
The simplest way of installing the current github version of the package
is opening an R session and executing the command:

     devtools::install_github(repo='amsstats', username='butwhywhy')

(This of course requires that you have the devtools package installed. 
We expect to be able to upload
the package to the official CRAN repository, so that the standard mechanisms
for installing packages in R work).

If you also want to play with the code, you can clone this git repository:

    git clone https://github.com/butwhywhy/amsstats.git
    cd amsstats

Modify whatever you want, and install executing

    devtools::install(pkg='.')

from within an R session in the amsstats home directory.

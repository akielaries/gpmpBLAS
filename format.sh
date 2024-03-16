#!/bin/sh
printf "[+] Formatting FORTRAN files...\n"
cd ../ && find . -type f -name "*.f90*" -exec fprettify {} --indent 4 \;

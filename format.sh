#!/bin/sh
printf "[+] Formatting FORTRAN files...\n"
cd ../ && find . -type f -name "*.f" -exec fprettify {} --indent 4 \;

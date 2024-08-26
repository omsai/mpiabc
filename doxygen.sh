#!/bin/sh
doxygen "$1" && sed -i '/\input{README_8md}/d' latex/refman.tex

#!/bin/sh

PYFILES="*.py biopepa/*.py facile/*.py BioModels_Database-r22-sbml_files/*.py web-biopepa-latex/*.py"

pylint -r no ${PYFILES}
while getopts ":p" opt; do
  case $opt in
    p)
      pychecker ${PYFILES}     
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done


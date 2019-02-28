## A script for adding this library to the Python path.

PYTHONEXE=python3

SITEDIR=`${PYTHONEXE} -c 'import site; site._script()' --user-site`
CURDIR=`pwd`

echo "Adding
   $CURDIR
to the python path in
   ${SITEDIR}/smilei.pth"

echo $CURDIR >> ${SITEDIR}/smilei.pth

version=$1

python setup.py sdist
twine upload dist/cfutils-${version}.tar.gz
python setup.py clean

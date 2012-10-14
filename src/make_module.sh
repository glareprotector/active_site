swig -c++ -python test.i 
g++ -fpic -c test.cpp test_wrap.cxx globals.cpp sample.cpp -I/usr/include/python2.6 -I /usr/local/lib/python2.6/dist-packages/mpi4py-1.3-py2.6-linux-x86_64.egg/mpi4py/include -g
g++ -shared globals.o test.o test_wrap.o /usr/lib/libpython2.6.so -o _test.so -g
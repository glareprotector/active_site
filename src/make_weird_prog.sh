
g++ -c test.cpp globals.cpp sample.cpp -I/usr/include/python2.6 -I /usr/local/lib/python2.6/dist-packages/mpi4py-1.3-py2.6-linux-x86_64.egg/mpi4py/include -g
g++ globals.o test.o /usr/lib/libpython2.6.so -o prog -g
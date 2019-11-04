This generates an executable compiled with g++ with
gperftools heap profiling enabled

to build enter
make
at the command prompt

gcc [...] -o myprogram -ltcmalloc
HEAPPROFILE=/tmp/profile ./myprogram

See

http://code.google.com/p/gperftools/wiki/GooglePerformanceTools

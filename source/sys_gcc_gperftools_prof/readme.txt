This generates an executable compiled with g++ with
gperftools profiling enabled

to build enter
make
at the command prompt

gcc [...] -o myprogram -lprofiler
CPUPROFILE=/tmp/profile ./myprogram

See

http://code.google.com/p/gperftools/wiki/GooglePerformanceTools

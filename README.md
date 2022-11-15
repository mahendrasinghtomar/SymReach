# SymReach
Additionally Requires:

    Gnuplot: http://gnuplot.info/
    Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
    Ginac: https://ginac.de/Download.html

Please edit makefile (in every example folder) to adjust directory path.

Libraries to be linked: -lboost_filesystem -lboost_iostreams -lboost_system -lginac -lcln

If compiler else than g++ then please edit ReachableSet2cpp.h: line 1654
        
       system("g++ -shared func_from_file.cpp -o func_from_file.so");
to adjust compiler.

Publication: https://easychair.org/publications/open/gjfh
    Fabian Immler, Matthias Althoff,[et al, including Mahendra Singh Tomar], ``ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics" 2018  ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems

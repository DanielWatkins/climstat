# climstat
Robust statistics for climate data

climstat is (to be) a Python library containing implementations of robust statistics useful for climate data analysis. The initial motivation was to have access to the statistical tools described in the paper "Resistant, robust, and Non-parametric Techniques for the Analysis of CLimate Data" by John R. Lazante (International Journal of CLimatology, 16, pp. 1197-1226, 1996).

climstat also includes tools for downloading and processing IGRA2 radiosonde data. (An updated version of this code was added in version 0.7 of Unidata's Siphon package as `IGRAUpperAir()`.)


climstat is for Python 3 and depends on Numpy and Scipy.stats
Currently, I am using Numpy 1.10.4, Scipy 0.18.1, and Python 3.4.5.

Yes, there is another package called climstats in Python and in R. No, I don't really care. If I get enough code put together on this that it could be published, I'll change the name.


# PyGBSThr 
---
``PyGBSTr`` is a python software that solves Gaussian Boson Sampling problem (threshold detection, without losses and noise) using the [A.S. Popova and A.N. Rubtsov method](https://arxiv.org/abs/2106.01445). 
``PyGBSTr`` contains both the exact calculation and the approximate approach.    

This project consists of two main parts: visual demonstration of our method on a small problem in ``demo.ipynb`` and the optimized realization of the approach for random matrices in  ``start_up.py``.

# Installation 
---
To install ``PyGBSTr`` you can download [Docker  Desktop](https://www.docker.com/get-started) to use our code more easily. 
If you would not like to use Docker, you also could follow the instructions below.

To run our ``demo.ipynb`` you need to install [Jupyter](https://jupyter.org/) enviroment, ``Python 3`` and the folloing libraries:
* numpy 
* matplotlib 
* random
* math
* itertools  
* [sympy](https://www.sympy.org/en/index.html)
* [thewalrus](https://the-walrus.readthedocs.io/en/latest/) 
* [strawberryfields](https://strawberryfields.ai/)

To run ``start_up.ipynb`` you need to install``g++`` compiler, ``Python 3`` and the folloing libraries:

* numpy 
* random
* math
* subprocess
* sys
 


# Usage 
---

We recommend you start with the ``demo.ipynb`` file to familiarize yourself with our method. It represents all meaningful steps of the computation using the small 8-mode GBS problem. To advanced use ``PyGBSTr`` the ``start_up.py`` file was made. 

### Docker
To run ``demo.ipynb`` you can use the following command
```
docker run -p 8080:8080 stacy8popova/py_gbs_th
```
To run ``start_up.ipynb`` you can insert
```
docker run -it -v $(pwd):/home/devel/in_out stacy8popova/py_gbs_thr python3 ./start_up.py
```
for Linux;

```
docker run -it -v ${pwd}:/home/devel/in_out stacy8popova/py_gbs_thr python3 ./start_up.py
```
for Windows PowerShell;

```
docker run -it -v %cd%:/home/devel/in_out stacy8popova/py_gbs_thr python3 ./start_up.py
```
for Windows command line (Cmd.exe).

### Linux

To run ``start_up.py`` you need to insert in a command line

```
make
python3 start_up.py
```
### Input

``start_up.py`` requires a user input 
```
Do you need the exact calculation? (yes/no): 
Number of modes:
Squeezing parameter:
Number of clicked detectors:
```
after that it launches ``.py`` and ``.cpp`` files sequentially.  

We recommend to start with ``Number of modes: m``, where ``m`` is 10 to 30, ``Squeezing parameter: r`` , where ``r`` is 1.4 to 1.8, ``Number of clicked detectors: n``, where ``n`` is 9 to 30, ``n/m``>0.5. 

The output is 

```
Finished: Matrix.py
Finished: Exact.cpp   (if the input 'yes')
Finished: Minors.cpp
Finished: Moments.py
Finished: Approximation.py
Finished: Result.py
```
and files in the directory ``\in_out``

* ``GBS_matrix.dat`` is symmetrical complex-valued matrix defining the multimode Gaussian state of the entire system;
* ``Initial_state.dat`` contains the parameters of the input single-mode vacuum states;  
* ``Parmeters_of_interferometer.dat`` contains the parameters of a random interferometer;
* ``Sample.dat`` contains the random string of clicked detectors numbers (measurement outcome);
* ``Submatrix.dat`` is a submatrix built from the rows and columns of ``Sample.dat``;
* ``Exact.dat`` is the column of the exact probabilities (unnormalized);
* ``Minors0-1.dat`` contains the partition functions for minors of the 0 and 1st order; 
* ``Minors2.dat`` contains the partition functions for minors of the 2nd order; 
* ``Minors3.dat`` contains the partition functions for minors of the 3rd order;
* ``Minors4.dat`` contains the partition functions for minors of the 4th order; 
* ``Moments.dat`` contains moments up to the 4th order;
* ``Approximation.dat`` contains approximate moments up to the 4th order;
* ``Result.dat`` contains the exact probability of the specific measurement outcome and estimated one for different orders of approximation.


# Support 
---

Please, don't hesitate to share your suggestions, comments or critics with the authors and send them to a.popova@rqc.ru.

### Warning
There are the range of parameters, where your code performs the best. At the same time, we think that this version can be improved to obtain stabler result of a wider range of parameters.
In ``demo.ipynb`` and ``start_up.py`` you can find some comments after ``#`` about meaning of a specific part of the code and ways to improve it. 

# Acknowledgments
---
We appreciate the [Xanadu](https://www.xanadu.ai/) team for the free and open-source libraries that enabled us to benchmark our method.

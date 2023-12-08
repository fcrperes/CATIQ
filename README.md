# CATIQ - CAT-state Injected Quantum Circuits

In this repository, you can find the Python code used to produce the results
presented in our paper "Non-stabilizerness and entanglement from cat-state
injection" [https://arxiv.org/abs/2305.19988].

Usage is fairly straighforward: Users need only open and make changes to the
file `main.py`.

Contact: filipa.peres@inl.int
-------------------------------------------------------------------------------


## Reproducing numerical experiments in the paper:

To reproduce the results in the paper, users should open the file `main.py` and
use the values presented below. The generation of each Figure should take
between 25min and 50min (depending on the experiment) in a standard laptop
computer. To reduce the run time, one can reduce the nr. of qubits, circuits,
or of evaluations of the OTOC.


* Common to all figures:
nq = 10
nlayers = 10
ncircs = 10 
evals = 100
rand = True
flagh = (either True or False: it is irrelevant!)
flags = (either True or False: it is irrelevant!)


* Figure 2(a):
gset = Clifford_separable
t_inject = []


* Figure 2(b):
gset = Clifford
t_inject = []


* Figure 2(c): ()
gset = Clifford_plus_T
t_inject = []


* Figure 5(a):
gset = Clifford_separable
t_inject = random.sample(range(5,41), 16)
t_inject.sort()


* Figure 5(b):
gset = Clifford_separable
t_inject = random.sample(range(5,91), 36)
t_inject.sort()



## Copyright:

Copyright (C) 2023  Filipa C. R. Peres

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

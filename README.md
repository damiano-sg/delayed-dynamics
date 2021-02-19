# delayed-dynamics

**Title of my Master thesis:** The effect of Delayed Dynamics on the Stability of Ecosystems. _http://tesi.cab.unipd.it/65108/_

**Abstract:**
Differential Equations and random matrix theory have found applications in many fields, ranging from physics, number theory and ecology. Fifty years ago, in a seminal article entitled "Will a Large Complex System be Stable?", Robert May showed that increasingly large or complex ecological networks have negligibly small probability to be stable. However, from field studies in ecology experimental evidences point on the opposite direction: large and complex ecosystems seem to be stable. This apparent contradiction is known as the Complexity-Stability paradox. May and other scientists analyzed large networks in which species interact at random and studied analytically the asymptotic stability of such systems as a function of the number of species and the number of interactions through the use of random matrix theory.
My master thesis is aimed to critically review these results and show how the introduction of delay in ecosystem population dynamics may have a strong impact on ecosystem dynamics and change important aspects of the Complexity-Stability paradox.


**Objectives:**
The main objective of this repository is to characterize the stability of the linearized system of Delay Differential Equations (DDEs) in the form:

x' (t) = A x(t) + B x(t-T)

We study how different delays change the distribution of the effective eigenvalues lambda of a generic DDE. In particular how the rightmost of them (the one which tell us if the system is stable or not) varies as a function of the parameters.

What we call eigevalues are the solutions of the characteristic equation obtained from the DDE by using the Linear Stability Analysis, i.e. by substituting solutions of the form:   x(t) = exp(\lambda t) v    in the DDE obtaining an equation for lambda.

Each of the following files is a jupyter notebook which displays the different characteristics which enter in the description of the stability of an ecosystem described by the linearized DDE above.

## Files
1. _RMG.py_ : class used to generate different types of random matrices.
2. _Distribution_eigenvalues.ipynb_ : plots of different types of eigenvalue's distributions.
3. _Analytical_results.ipynb_ : analytical study of the eigenvalue equation and how that is connected to the delay term.
4. _Maximal_eigenvalue.ipynb_ : scaling laws of the real part of the rightmost eigenvalue as a function of delay and other parameters.
5. _Stability_trajectories.ipynb_ : comparison between the stability relations obtained using the rightmost eigenvalue technique and the actual stability of the solutions obtained by observing the convergence to equilibrium of the trajectories.
6. _Other_branches.ipynb_ : check that the other branches of the Lambert function don't influence the rightmost eigenvalue (which remains always in the principal branch)
7. _Contour_distributions.ipynb_ : check that the borders of the distributions are mapped exactly from the borders of the Random matrices eigenvalue's distributions.
8. _Distributed_delay.ipynb_ : same analysis of the previous cases but with a different type of delay entering the DDE: the delay is now distributed, i.e. now the DDE is an integro-differential equation. The main results are the same as the fixed-delay case

